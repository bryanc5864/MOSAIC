"""Post-hoc crop the aligned-latent UMAPs to zoom into the data region.

The original matplotlib axes span all data including outliers, squishing the
main cluster blob into a corner of panels 2 and 3. We can't regenerate
without per-cell latent .npy files (not committed). As a readability fix we
crop each panel to its data-rich region and reassemble side-by-side.
"""
from pathlib import Path
import numpy as np
from PIL import Image

FIG = Path(__file__).resolve().parent.parent / "figures"


def crop_scatter_dots(arr, pad=5):
    """Return bbox of colored (non-white, non-gray) pixels — i.e., scatter dots.

    Axis lines, ticks, text, and titles are grayscale and should be excluded.
    """
    r, g, b = arr[..., 0].astype(int), arr[..., 1].astype(int), arr[..., 2].astype(int)
    max_c = np.maximum.reduce([r, g, b])
    min_c = np.minimum.reduce([r, g, b])
    saturation = max_c - min_c  # simple chroma proxy
    not_white = max_c < 245
    colored = (saturation > 25) & not_white
    if not colored.any():
        return None
    ys, xs = np.where(colored)
    # Use 2nd–98th percentile to ignore stray outlier dots pulling the bbox wide
    x0 = int(np.percentile(xs, 2)); x1 = int(np.percentile(xs, 98))
    y0 = int(np.percentile(ys, 2)); y1 = int(np.percentile(ys, 98))
    return (max(x0 - pad, 0), max(y0 - pad, 0),
            min(x1 + pad, arr.shape[1] - 1), min(y1 + pad, arr.shape[0] - 1))


def rebuild_umap_fig(in_path: Path, out_path: Path):
    img = Image.open(in_path).convert("RGB")
    W, H = img.size
    arr = np.asarray(img)
    print(f"  input {W}x{H}")

    # The image has: suptitle at top, 3 panels side by side, legend at right.
    # Layout from fig, axes = plt.subplots(1, 3, figsize=(14, 4.5)) then legend
    # anchored at (1.03, 0.5). Approximate horizontal splits:
    # title band  ~  top 8%
    # panels      ~  (title_band) .. (just before legend ~ 85% of W)
    # legend      ~  right ~15% of W
    title_h = int(H * 0.08)
    legend_x = int(W * 0.86)

    # Isolate panels strip
    panels_strip = arr[title_h:, :legend_x, :]
    legend_strip = arr[:, legend_x:, :]
    title_strip = arr[:title_h, :legend_x, :]
    panel_w = panels_strip.shape[1] // 3
    ps_h = panels_strip.shape[0]

    # Collect per-panel bboxes on scatter dots only
    bboxes = []
    raw_panels = []
    for i in range(3):
        panel = panels_strip[:, i * panel_w:(i + 1) * panel_w, :]
        raw_panels.append(panel)
        bboxes.append(crop_scatter_dots(panel, pad=30))

    # Use the UNION of y-bounds and per-panel x-bounds so all panels share a
    # common height (keeps panel titles and axes consistent), while x is
    # panel-specific so each panel zooms into its own horizontal data spread.
    valid_bboxes = [b for b in bboxes if b is not None]
    if valid_bboxes:
        y0_glob = min(b[1] for b in valid_bboxes)
        y1_glob = max(b[3] for b in valid_bboxes)
    else:
        y0_glob, y1_glob = 0, panels_strip.shape[0] - 1

    cropped_panels = []
    for i, (panel, bbox) in enumerate(zip(raw_panels, bboxes)):
        if bbox is None:
            x0, x1 = 0, panel.shape[1] - 1
        else:
            x0, x1 = bbox[0], bbox[2]
        cp = panel[y0_glob:y1_glob + 1, x0:x1 + 1, :]
        cropped_panels.append(Image.fromarray(cp))
        print(f"  panel {i} bbox {bbox} -> {cp.shape[1]}x{cp.shape[0]}")

    # Fit every panel into a uniform box preserving aspect ratio (no squashing);
    # pad remaining area with white so panels stay visually consistent.
    target_h = max(p.size[1] for p in cropped_panels)
    target_w = max(p.size[0] for p in cropped_panels)
    resized = []
    for p in cropped_panels:
        scale = min(target_w / p.size[0], target_h / p.size[1])
        new_w = int(p.size[0] * scale); new_h = int(p.size[1] * scale)
        scaled = p.resize((new_w, new_h), Image.LANCZOS)
        canvas = Image.new("RGB", (target_w, target_h), "white")
        canvas.paste(scaled, ((target_w - new_w) // 2, (target_h - new_h) // 2))
        resized.append(canvas)

    # Padding between panels
    gap = 40
    panels_w = sum(p.size[0] for p in resized) + gap * (len(resized) - 1)
    panels_img = Image.new("RGB", (panels_w, target_h), "white")
    x = 0
    for p in resized:
        panels_img.paste(p, (x, 0))
        x += p.size[0] + gap

    # Scale legend to panel height
    legend_im = Image.fromarray(legend_strip)
    lscale = target_h / legend_im.size[1]
    legend_im = legend_im.resize((int(legend_im.size[0] * lscale), target_h), Image.LANCZOS)

    # Scale title to new panels width
    title_im = Image.fromarray(title_strip)
    tscale = panels_w / title_im.size[0]
    title_im = title_im.resize((panels_w, int(title_im.size[1] * tscale)), Image.LANCZOS)

    out_w = panels_w + legend_im.size[0]
    out_h = title_im.size[1] + target_h
    out = Image.new("RGB", (out_w, out_h), "white")
    out.paste(title_im, (0, 0))
    out.paste(panels_img, (0, title_im.size[1]))
    out.paste(legend_im, (panels_w, title_im.size[1]))
    out.save(out_path, optimize=True)
    print(f"  saved {out_path.name} {out.size}")


if __name__ == "__main__":
    for name in ("fig1_aligned_latent_pbmc10k_multiome",
                 "fig1_aligned_latent_brain3k_multiome"):
        print(name)
        rebuild_umap_fig(FIG / f"{name}.png", FIG / f"{name}.png")
