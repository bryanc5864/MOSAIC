# MIT License
# Part of MOSAIC
"""Preprocessing for the paired multi-omics benchmarks.

Load the 10x h5, split by modality, QC jointly so the pairing survives, then
RNA: normalize_total -> log1p -> 2000 HVGs; ATAC: binarize -> TF-IDF -> 50 LSI
dims with the first dropped. Cross-modal targets go under .obsm so the training
loop never has to touch the partner modality's full matrix. Output lands in
data/processed/<dataset_id>_{rna,atac}.h5ad.

Multiome and CITE-seq have separate loaders.
"""

from __future__ import annotations

import gc
import hashlib
import warnings
from dataclasses import asdict
from pathlib import Path

import anndata as ad
import numpy as np
import scanpy as sc
import scipy.sparse as sp
from sklearn.decomposition import TruncatedSVD

from src.data.datasets import DATASETS, DatasetSpec
from src.utils.paths import PROCESSED_DIR

warnings.filterwarnings("ignore", category=UserWarning, module="anndata")
warnings.filterwarnings("ignore", category=FutureWarning, module="scanpy")


RNA_MIN_GENES = 200              # cells must express at least this many genes
RNA_MIN_CELLS_PER_GENE = 3       # genes expressed in at least this many cells
RNA_N_HVG = 2000                 # number of highly variable genes to keep
RNA_NORM_TARGET = 1e4            # total-count normalization target
RNA_PCA_DIM = 50                 # RNA PCA dim used for clustering + cross-modal target

ATAC_MIN_PEAKS = 1000            # cells must have at least this many accessible peaks
ATAC_MIN_CELLS_PER_PEAK = 10     # peaks accessible in at least this many cells
# 10K rather than 50K peaks: keeps the decoder's output layer tractable on one
# RTX 3080 and still covers the variable chromatin landscape.
ATAC_N_VAR_PEAKS = 10_000
ATAC_LSI_DIM = 50                # total LSI dims kept (first one dropped downstream)

LEIDEN_RES = 0.8                 # leiden resolution for cluster labeling

SEED = 0


# tf-idf + lsi inline rather than episcanpy, which kept breaking the env


def tfidf(X: sp.spmatrix) -> sp.spmatrix:
    """TF-IDF on a sparse cells x peaks matrix, Signac/ArchR formulation.

    Input is binarized first, so this is presence/absence of accessibility.
    """
    X = (X > 0).astype(np.float32)           # binarize
    X = sp.csr_matrix(X)
    n_cells = X.shape[0]
    row_sums = np.asarray(X.sum(axis=1)).ravel()
    row_sums[row_sums == 0] = 1.0            # avoid div by zero
    tf_inv = 1.0 / row_sums
    X_tf = X.multiply(tf_inv[:, None]).tocsr()
    peak_sums = np.asarray(X.sum(axis=0)).ravel()
    idf = np.log(1.0 + (n_cells / np.maximum(peak_sums, 1.0)))
    X_tfidf = X_tf.multiply(idf[None, :]).tocsr()
    return X_tfidf


def lsi(X_tfidf: sp.spmatrix, n_components: int, seed: int = SEED) -> np.ndarray:
    """Truncated SVD on the TF-IDF matrix.

    Caller should drop component 1, which tracks read depth.
    """
    svd = TruncatedSVD(n_components=n_components, random_state=seed, algorithm="arpack")
    Z = svd.fit_transform(X_tfidf)
    return Z.astype(np.float32)


def _md5_of_array(X) -> str:
    """md5 of the array bytes, for integrity logging"""
    if sp.issparse(X):
        data = np.ascontiguousarray(X.data).tobytes() + \
               np.ascontiguousarray(X.indices).tobytes() + \
               np.ascontiguousarray(X.indptr).tobytes()
    else:
        data = np.ascontiguousarray(X).tobytes()
    return hashlib.md5(data).hexdigest()


def load_10x_multiome_h5(path: Path) -> tuple[ad.AnnData, ad.AnnData]:
    """split a 10x multiome h5 into (rna, atac) with matching row order"""
    full = sc.read_10x_h5(path, gex_only=False)
    full.var_names_make_unique()
    ft = full.var["feature_types"].astype(str).values
    rna_mask = ft == "Gene Expression"
    atac_mask = ft == "Peaks"
    if rna_mask.sum() == 0 or atac_mask.sum() == 0:
        raise ValueError(
            f"Expected both 'Gene Expression' and 'Peaks' feature types; "
            f"got {np.unique(ft)}"
        )
    rna = full[:, rna_mask].copy()
    atac = full[:, atac_mask].copy()
    assert np.array_equal(rna.obs_names.values, atac.obs_names.values), \
        "RNA and ATAC cell orders must match after split"
    # survives later subsetting; this is the ground-truth pairing for eval
    rna.obs["pair_idx"] = np.arange(rna.n_obs)
    atac.obs["pair_idx"] = np.arange(atac.n_obs)
    return rna, atac


def qc_filter_joint(rna: ad.AnnData, atac: ad.AnnData) -> tuple[ad.AnnData, ad.AnnData]:
    """drop cells failing either modality's QC, keeping the pairing intact"""
    n0 = rna.n_obs
    rna_n_genes = np.asarray((rna.X > 0).sum(axis=1)).ravel()
    atac_n_peaks = np.asarray((atac.X > 0).sum(axis=1)).ravel()
    keep = (rna_n_genes >= RNA_MIN_GENES) & (atac_n_peaks >= ATAC_MIN_PEAKS)
    print(f"  [qc] cells: {n0} -> {keep.sum()}  "
          f"(rna_min_genes>={RNA_MIN_GENES}, atac_min_peaks>={ATAC_MIN_PEAKS})")
    rna_f = rna[keep].copy()
    atac_f = atac[keep].copy()
    n_genes_before = rna_f.n_vars
    sc.pp.filter_genes(rna_f, min_cells=RNA_MIN_CELLS_PER_GENE)
    n_peaks_before = atac_f.n_vars
    peak_prevalence = np.asarray((atac_f.X > 0).sum(axis=0)).ravel()
    atac_f = atac_f[:, peak_prevalence >= ATAC_MIN_CELLS_PER_PEAK].copy()
    print(f"  [qc] genes: {n_genes_before} -> {rna_f.n_vars}, "
          f"peaks: {n_peaks_before} -> {atac_f.n_vars}")
    return rna_f, atac_f


def preprocess_rna(rna: ad.AnnData) -> ad.AnnData:
    """normalize -> log1p -> hvg -> pca"""
    # zinb decoder needs the raw counts
    rna.layers["counts"] = rna.X.copy()
    sc.pp.normalize_total(rna, target_sum=RNA_NORM_TARGET)
    sc.pp.log1p(rna)
    sc.pp.highly_variable_genes(rna, n_top_genes=RNA_N_HVG, flavor="seurat_v3",
                                 layer="counts", subset=True)
    sc.pp.scale(rna, max_value=10)
    sc.tl.pca(rna, n_comps=RNA_PCA_DIM, random_state=SEED)
    # leiden gives us the cell-type labels
    sc.pp.neighbors(rna, n_neighbors=15, n_pcs=RNA_PCA_DIM, random_state=SEED)
    sc.tl.leiden(rna, resolution=LEIDEN_RES, random_state=SEED, flavor="igraph", n_iterations=2)
    sc.tl.umap(rna, random_state=SEED)
    return rna


def preprocess_atac(atac: ad.AnnData) -> ad.AnnData:
    """variable peaks -> binarize -> tf-idf -> lsi"""
    # bce decoder reads this layer
    atac.layers["binary"] = (atac.X > 0).astype(np.float32).tocsr() \
        if sp.issparse(atac.X) else (atac.X > 0).astype(np.float32)
    prev = np.asarray((atac.X > 0).sum(axis=0)).ravel().astype(np.float32)
    # bernoulli variance as the ranking proxy
    n_cells = atac.n_obs
    p_hat = prev / n_cells
    var_proxy = p_hat * (1 - p_hat)
    n_keep = min(ATAC_N_VAR_PEAKS, atac.n_vars)
    top_peaks = np.argsort(-var_proxy)[:n_keep]
    atac = atac[:, np.sort(top_peaks)].copy()
    print(f"  [atac] kept top {atac.n_vars} variable peaks")
    X_tfidf = tfidf(atac.X)
    Z_lsi = lsi(X_tfidf, n_components=ATAC_LSI_DIM + 1, seed=SEED)
    # drop component 1, it tracks depth
    Z_lsi_final = Z_lsi[:, 1:1 + ATAC_LSI_DIM].astype(np.float32)
    atac.obsm["X_lsi"] = Z_lsi_final
    # encoder input is log1p of tf-idf, straight off .X
    atac.layers["tfidf"] = X_tfidf.tocsr()
    atac.X = sp.csr_matrix(np.log1p(X_tfidf.toarray())).astype(np.float32)
    atac.obsm["X_pca"] = Z_lsi_final  # alias so scanpy's neighbors/umap work
    sc.pp.neighbors(atac, n_neighbors=15, use_rep="X_lsi", random_state=SEED)
    sc.tl.umap(atac, random_state=SEED)
    return atac


def add_cross_modal_targets(rna: ad.AnnData, atac: ad.AnnData) -> None:
    """Write cross-modal prediction targets into each modality's .obsm.

    Targets are cluster centroids, not per-cell partner features. For a cell in
    modality A, y_cross is the mean of the partner modality's low-dim summary
    over all cells in the same leiden cluster.

    Per-cell targets were the original design and they failed: the head memorized
    10K unique fingerprints, train_pred went to 0 while val_pred climbed (run-002).
    Sharing one target across a cluster forces it to learn cluster identity, which
    generalizes to held-out cells and is still genuinely cross-modal.

    Caveat worth knowing: the clusters come from leiden on RNA and are copied onto
    ATAC using the paired ground truth. Fine on a paired benchmark, but a truly
    unpaired run would have to cluster one modality alone or bring external labels.

    Centroids are standardized per dimension over cells, not over centroids, so
    the MSE stays on the old scale.
    """
    assert rna.n_obs == atac.n_obs, "Must be called on paired, QC-matched AnnDatas"
    # copy leiden labels from RNA onto ATAC
    atac.obs["leiden"] = rna.obs["leiden"].values
    atac.obs["cell_type"] = rna.obs["leiden"].values
    rna.obs["cell_type"] = rna.obs["leiden"].values

    rna_pca = np.asarray(rna.obsm["X_pca"], dtype=np.float32)
    atac_lsi = np.asarray(atac.obsm["X_lsi"], dtype=np.float32)
    rna_pca_std = (rna_pca - rna_pca.mean(0)) / (rna_pca.std(0) + 1e-6)
    atac_lsi_std = (atac_lsi - atac_lsi.mean(0)) / (atac_lsi.std(0) + 1e-6)

    clusters = rna.obs["leiden"].astype(str).values
    unique_clusters = np.unique(clusters)
    rna_centroid = {}
    atac_centroid = {}
    for c in unique_clusters:
        mask = clusters == c
        rna_centroid[c] = rna_pca_std[mask].mean(axis=0)
        atac_centroid[c] = atac_lsi_std[mask].mean(axis=0)

    # each cell gets its cluster's centroid in the other modality's space
    rna_y_cross = np.stack([atac_centroid[c] for c in clusters]).astype(np.float32)
    atac_y_cross = np.stack([rna_centroid[c] for c in clusters]).astype(np.float32)
    rna.obsm["y_cross"] = rna_y_cross
    atac.obsm["y_cross"] = atac_y_cross
    # kept for the entropy correlation diagnostics
    rna.obsm["X_pca_std"] = rna_pca_std
    atac.obsm["X_lsi_std"] = atac_lsi_std


def save_processed(dataset_id: str, rna: ad.AnnData, atac: ad.AnnData) -> dict:
    """write both AnnDatas, return integrity metadata"""
    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
    rna_path = PROCESSED_DIR / f"{dataset_id}_rna.h5ad"
    atac_path = PROCESSED_DIR / f"{dataset_id}_atac.h5ad"
    # anndata chokes on the sparse tfidf layer sometimes
    if "tfidf" in atac.layers:
        del atac.layers["tfidf"]
    rna.write_h5ad(rna_path, compression="gzip")
    atac.write_h5ad(atac_path, compression="gzip")
    meta = {
        "dataset_id": dataset_id,
        "rna_path": str(rna_path),
        "atac_path": str(atac_path),
        "rna_shape": list(rna.shape),
        "atac_shape": list(atac.shape),
        "rna_md5": _md5_of_array(rna.X),
        "atac_md5": _md5_of_array(atac.X),
        "n_clusters": int(rna.obs["leiden"].nunique()),
        "cluster_sizes": rna.obs["leiden"].value_counts().to_dict(),
    }
    return meta


def load_10x_citeseq_h5(path: Path) -> tuple[ad.AnnData, ad.AnnData]:
    """split a 10x cite-seq h5 into (rna, protein)"""
    full = sc.read_10x_h5(path, gex_only=False)
    full.var_names_make_unique()
    ft = full.var["feature_types"].astype(str).values
    rna_mask = ft == "Gene Expression"
    prot_mask = ft == "Antibody Capture"
    if rna_mask.sum() == 0 or prot_mask.sum() == 0:
        raise ValueError(
            f"Expected 'Gene Expression' and 'Antibody Capture'; got {np.unique(ft)}"
        )
    rna = full[:, rna_mask].copy()
    prot = full[:, prot_mask].copy()
    assert np.array_equal(rna.obs_names.values, prot.obs_names.values), \
        "RNA and protein cell orders must match"
    rna.obs["pair_idx"] = np.arange(rna.n_obs)
    prot.obs["pair_idx"] = np.arange(prot.n_obs)
    return rna, prot


def preprocess_protein(prot: ad.AnnData, drop_isotype: bool = True) -> ad.AnnData:
    """CLR-normalize CITE-seq ADT counts.

    Protein counts are dense, so CLR rather than TF-IDF. Isotype controls carry
    no cell-type signal and are dropped by default. PCA is only 16-dim because
    there are only ~14 real markers.
    """
    if drop_isotype:
        is_control = np.array(["control" in str(n).lower() or "IgG" in str(n)
                                for n in prot.var_names])
        print(f"  [prot] dropping {is_control.sum()} isotype controls; "
              f"keeping {(~is_control).sum()} protein markers")
        prot = prot[:, ~is_control].copy()

    prot.layers["counts"] = prot.X.copy()
    # clr: log(x / geometric mean(x)) per cell
    X = np.asarray(prot.X.toarray() if sp.issparse(prot.X) else prot.X,
                   dtype=np.float64)
    X_p = X + 1.0
    log_X = np.log(X_p)
    geo_mean = log_X.mean(axis=1, keepdims=True)
    clr = log_X - geo_mean
    prot.X = clr.astype(np.float32)
    n_comps = min(16, prot.n_vars - 1)
    sc.tl.pca(prot, n_comps=n_comps, random_state=SEED)
    # alias so the cross-modal-target code doesn't need a protein special case
    prot.obsm["X_lsi"] = prot.obsm["X_pca"]
    sc.pp.neighbors(prot, n_neighbors=15, use_rep="X_pca", random_state=SEED)
    sc.tl.umap(prot, random_state=SEED)
    return prot


def run_pipeline_citeseq(dataset_id: str) -> dict:
    """end-to-end cite-seq preprocessing"""
    spec = DATASETS[dataset_id]
    print(f"\n[preprocess] {dataset_id}")
    print(f"  source: {spec.raw_path}")
    rna, prot = load_10x_citeseq_h5(spec.raw_path)
    print(f"  loaded: rna {rna.shape}, protein {prot.shape}")

    # no protein threshold: every cell expresses some of every marker
    n0 = rna.n_obs
    rna_n_genes = np.asarray((rna.X > 0).sum(axis=1)).ravel()
    keep = rna_n_genes >= RNA_MIN_GENES
    print(f"  [qc] cells: {n0} -> {keep.sum()}")
    rna = rna[keep].copy()
    prot = prot[keep].copy()
    sc.pp.filter_genes(rna, min_cells=RNA_MIN_CELLS_PER_GENE)

    print("  [rna] running scanpy pipeline")
    rna = preprocess_rna(rna)
    print("  [prot] running CLR + PCA")
    prot = preprocess_protein(prot)

    # works unchanged thanks to the X_lsi alias set above
    add_cross_modal_targets(rna, prot)
    # protein is saved under the "atac" suffix so training can keep using
    # recon_layer='binary'; binarize each marker against the per-cell median
    prot_dense = np.asarray(prot.X)
    median_per_cell = np.median(prot_dense, axis=1, keepdims=True)
    prot.layers["binary"] = (prot_dense > median_per_cell).astype(np.float32)

    meta = save_processed(dataset_id, rna, prot)
    print(f"  [save] {meta['rna_path']}")
    print(f"  [save] {meta['atac_path']}")
    print(f"  [save] {meta['n_clusters']} leiden clusters")
    gc.collect()
    return meta


def run_pipeline_10x_multiome(dataset_id: str) -> dict:
    """end-to-end multiome preprocessing"""
    spec = DATASETS[dataset_id]
    print(f"\n[preprocess] {dataset_id}")
    print(f"  source: {spec.raw_path}")
    rna, atac = load_10x_multiome_h5(spec.raw_path)
    print(f"  loaded: rna {rna.shape}, atac {atac.shape}")
    rna, atac = qc_filter_joint(rna, atac)
    print("  [rna] running scanpy pipeline")
    rna = preprocess_rna(rna)
    print("  [atac] running TF-IDF + LSI")
    atac = preprocess_atac(atac)
    add_cross_modal_targets(rna, atac)
    meta = save_processed(dataset_id, rna, atac)
    print(f"  [save] {meta['rna_path']}")
    print(f"  [save] {meta['atac_path']}")
    print(f"  [save] {meta['n_clusters']} leiden clusters")
    gc.collect()
    return meta


if __name__ == "__main__":
    import argparse
    import json

    p = argparse.ArgumentParser(description="preprocess a benchmark dataset")
    p.add_argument("dataset", help="Dataset id (e.g. pbmc10k_multiome)")
    args = p.parse_args()

    if args.dataset not in DATASETS:
        raise SystemExit(f"unknown dataset: {args.dataset}")

    spec = DATASETS[args.dataset]
    if "protein" in spec.modalities:
        meta = run_pipeline_citeseq(args.dataset)
    else:
        meta = run_pipeline_10x_multiome(args.dataset)
    meta_path = PROCESSED_DIR / f"{args.dataset}_meta.json"
    with meta_path.open("w") as f:
        json.dump(meta, f, indent=2, default=str)
    print(f"\n[done] metadata saved to {meta_path}")
