"""Rendering and previews for a detection result.

Everything here works on the smallest pyramid level: the overview image and the
per-box crops are drawn and cut from the same thumbnail, so box coordinates and
tissue masks line up without any resizing.

``preview_detection`` returns the matplotlib figures it builds, so it can be
used directly in a notebook without writing anything to disk. The matplotlib
backend is intentionally *not* forced here -- under Jupyter's inline backend the
returned figures display automatically. The command-line entry point selects the
non-interactive Agg backend itself for headless/cluster use.
"""

import os

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np
import skimage.measure

from .geometry import _other_tissue_mask
from .pyramid import pyramid_levels

_COMPOSITE_COLORS = np.array([
    [0.0, 0.6, 1.0],
    [1.0, 0.5, 0.0],
    [0.0, 1.0, 0.3],
    [1.0, 0.2, 0.5],
    [0.9, 0.9, 0.0],
    [0.6, 0.4, 1.0],
])


def _stretch(array, p_low=1.0, p_high=99.0):
    array = np.asarray(array, dtype=float)
    lo, hi = np.percentile(array, (p_low, p_high))
    if hi <= lo:
        lo, hi = float(array.min()), float(array.max())
    return np.clip((array - lo) / (hi - lo + 1e-9), 0.0, 1.0)


def _to_rgb(image):
    """Convert ``(c, y, x)`` or ``(y, x)`` data to display-ready RGB."""
    image = np.asarray(image)
    if image.ndim == 2:
        gray = _stretch(image)
        return np.dstack([gray, gray, gray])
    if image.ndim != 3:
        raise ValueError(f"Expected (c, y, x) or (y, x), got shape {image.shape}")

    rgb = np.zeros((*image.shape[1:], 3), dtype=float)
    for channel, data in enumerate(image):
        rgb += _stretch(data)[..., None] * _COMPOSITE_COLORS[
            channel % len(_COMPOSITE_COLORS)
        ]
    return np.clip(rgb, 0.0, 1.0)


def _tight_figure(rgb, dpi):
    """A figure whose single axes fills the canvas, one output pixel per input pixel."""
    height, width = rgb.shape[:2]
    fig = plt.figure(figsize=(width / dpi, height / dpi), dpi=dpi)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.imshow(rgb)
    ax.set_axis_off()
    return fig, ax


def _draw_mask_border(ax, mask, offset_yx, color, linewidth):
    y0, x0 = offset_yx
    contours = skimage.measure.find_contours(np.pad(mask, 1).astype(float), 0.5)
    for contour in contours:
        contour = contour - 1
        ax.plot(
            contour[:, 1] + x0,
            contour[:, 0] + y0,
            color=color,
            linewidth=linewidth,
        )


def overview_figure(
    image,
    boxes,
    dpi=150,
    included_color="lime",
    subtract_color="red",
    box_color="cyan",
    linewidth=1.3,
    fontsize=11,
):
    """Build and return the overview figure: included/subtracted tissue and boxes.

    ``image`` is the smallest pyramid level as an array. The figure is returned
    open (not saved, not closed) so the caller can display or save it.
    """
    rgb = _to_rgb(image)
    fig, ax = _tight_figure(rgb, dpi)

    for i, box in enumerate(boxes, start=1):
        y0, x0, y1, x1 = box.bbox_thumb
        ax.add_patch(
            Rectangle(
                (x0, y0),
                x1 - x0,
                y1 - y0,
                fill=False,
                edgecolor=box_color,
                linewidth=linewidth,
            )
        )

        _draw_mask_border(
            ax,
            box.mask,
            offset_yx=(y0, x0),
            color=included_color,
            linewidth=linewidth,
        )

        others = _other_tissue_mask(box, boxes)
        if others.any():
            _draw_mask_border(
                ax,
                others,
                offset_yx=(y0, x0),
                color=subtract_color,
                linewidth=linewidth,
            )

        ax.text(
            x0 + 2,
            y0 + 2,
            str(i),
            color="black",
            fontsize=fontsize,
            va="top",
            ha="left",
            fontweight="bold",
            bbox=dict(facecolor=box_color, edgecolor="none", pad=1.5, alpha=0.9),
        )

    return fig


def crop_figure(sdata, box, boxes, image_key="morphology_focus", dpi=150):
    """Build and return the crop figure for a single ``box``.

    The full bounding box is kept and only tissue assigned to other detected
    boxes is blacked out. Cut from the smallest pyramid level.
    """
    smallest = pyramid_levels(sdata, image_key)[0]
    y0, x0, y1, x1 = box.bbox_thumb
    crop = np.asarray(smallest.isel(y=slice(y0, y1), x=slice(x0, x1)))
    rgb = _to_rgb(crop)
    # ``_other_tissue_mask`` is already at smallest-level resolution and matches
    # the crop shape exactly, so no resizing is needed.
    rgb[_other_tissue_mask(box, boxes)] = 0.0
    fig, _ = _tight_figure(rgb, dpi)
    return fig


def preview_detection(
    sdata,
    boxes,
    outdir=None,
    image_key="morphology_focus",
    overview_name="tissue_detection.png",
    crop_prefix="tissue",
    crop_suffix="_preview",
    dpi=150,
    included_color="lime",
    subtract_color="red",
    box_color="cyan",
    save_overview=True,
    save_crops=True,
):
    """Build detection-preview figures and, if ``outdir`` is given, save them.

    ``boxes`` should come from ``detect_tissue(sdata)``. Both the overview and
    the crops are drawn from the smallest pyramid level.

    The overview (built when ``save_overview`` is true) shows:
      - ``included_color``: tissue belonging to each detected box
      - ``subtract_color``: neighbouring tissue removed from that box
      - ``box_color``: bounding boxes

    The crops (built when ``save_crops`` is true) keep the full rectangle and
    black out only tissue assigned to other detected boxes.

    Saving to disk is optional: pass ``outdir`` to write PNGs, or leave it as
    ``None`` to build the figures in memory only -- handy in a notebook, where
    the returned figures display inline. The figures are returned open (not
    closed); close them with ``matplotlib.pyplot.close`` when processing many
    images in a loop.

    Returns ``(overview_fig, crop_figs)``. ``overview_fig`` is ``None`` when the
    overview is disabled, and ``crop_figs`` is an empty list when crops are
    disabled.
    """
    smallest = pyramid_levels(sdata, image_key)[0]
    if outdir is not None:
        os.makedirs(outdir, exist_ok=True)

    overview_fig = None
    if save_overview:
        overview_fig = overview_figure(
            np.asarray(smallest),
            boxes,
            dpi=dpi,
            included_color=included_color,
            subtract_color=subtract_color,
            box_color=box_color,
        )
        if outdir is not None:
            overview_fig.savefig(os.path.join(outdir, overview_name), dpi=dpi)

    crop_figs = []
    if save_crops:
        for i, box in enumerate(boxes, start=1):
            fig = crop_figure(sdata, box, boxes, image_key=image_key, dpi=dpi)
            crop_figs.append(fig)
            if outdir is not None:
                fig.savefig(
                    os.path.join(outdir, f"{crop_prefix}_{i}{crop_suffix}.png"),
                    dpi=dpi,
                )

    if outdir is not None:
        wrote = []
        if save_overview:
            wrote.append("overview")
        if save_crops:
            wrote.append(f"{len(crop_figs)} crop(s)")
        print(f"Wrote {' + '.join(wrote) or 'nothing'} to {outdir}")

    return overview_fig, crop_figs
