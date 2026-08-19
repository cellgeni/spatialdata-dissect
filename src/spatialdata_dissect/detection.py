"""Detect the main tissue pieces on the smallest pyramid level.

This module also defines the two data structures the rest of the package is
built around: :class:`TissuePolicy` (detection parameters) and
:class:`TissueBox` (a detected component and its bounding box).

The public entry points are :func:`detect_tissue`, which returns largest-first
``TissueBox`` objects for a SpatialData image, and :func:`merge_boxes`, which
fuses several detected boxes into one.
"""

from dataclasses import dataclass

import numpy as np
import shapely
from scipy import ndimage

# Import submodules explicitly so they are available regardless of whether the
# installed scikit-image uses lazy submodule loading.
import skimage.filters
import skimage.measure
import skimage.morphology

from .pyramid import pyramid_levels


@dataclass(frozen=True)
class TissuePolicy:
    """Parameters used to detect the main tissue pieces on the smallest image level."""

    # it's the minimum component area for something to count as a real tissue piece
    # the same value also sets the hole-filling size in the support step 
    # raise it to ignore small fragments and fill larger interior gaps; lower it to keep tiny pieces and preserve small holes
    min_area_fraction: float = 0.002
    # morphological closing (dilate then erode) with a disk of this radius on the seed
    # it fills pin-holes and bridges bright pixels that are almost touching, so a speckled dense region becomes one solid blob
    # raise it to consolidate a broken-up seed; lower it to keep fine structure separate
    close_radius_px: int = 4
    # runs an opening (erode then dilate) right after it deletes isolated specks smaller than the disk
    # raise it if picking up scattered noise; set it to 0 to keep every bright speck
    open_radius_px: int = 2
    # gaussian-blurs the grayscale before thresholding it's also the sigma used by the density path below
    # larger values recover sparser tissue but smear intensities down; 0 thresholds the raw image
    support_blur_sigma_px: float = 4.0
    # sets the support threshold to at least this fraction of the Otsu seed threshold
    # the primary "how faint can included tissue be" 
    # lower includes fainter tissue; too low starts pulling in background
    support_threshold_fraction: float = 0.15
    # it stops the support threshold from sinking into background noise
    # on clean black backgrounds it's usually inert; on noisy backgrounds it's what prevents flooding
    # raise it to be more conservative near noise, lower it to allow the threshold closer to background
    support_bg_sigmas: float = 3.0
    # when positive it counts the local fraction of above-threshold pixels and includes any neighborhood
    # where that fraction exceeds this value even if no single pixel is bright. it's a fraction between 0 and 1
    support_density_threshold: float = 0.05
    # this join pieces that are close by closing the combined support mask with a disk of this radius
    # larger bridges more distant fragments; too large can fuse things that should stay separate
    support_close_radius_px: int = 8
    # pads each piece's bounding box by this fraction of its own width and height on every side
    # it only affects the crop rectangle, not the tissue mask
    # raise it if crops are clipping tissue at the edges
    box_margin_fraction: float = 0.05
    # grow each piece's box to swallow detached fragments within this many
    # smallest-level px; 0 = off. Only claims strays (non-kept tissue)
    # so it never merges two separate detected pieces.
    attach_stray_radius_px: int = 25
    # caps how many pieces are returned, after sorting largest-first
    # if you have more physical sections than this you'll lose the smallest ones
    # if you're getting spurious extra then lowering it trims to the biggest few
    max_candidates: int = 12


@dataclass(frozen=True)
class TissueBox:
    """A detected tissue component and its bounding box.

    ``bbox_thumb`` and ``mask`` live on the smallest pyramid level.
    ``scale_yx`` maps those coordinates to level-0 pixel coordinates.
    """

    bbox_thumb: tuple[int, int, int, int]  # y0, x0, y1, x1
    mask: np.ndarray
    scale_yx: tuple[float, float]
    area_px: int

    def rect_polygon(self):
        """Bounding-box rectangle in level-0 pixel coordinates."""
        y0, x0, y1, x1 = self.bbox_thumb
        sy, sx = self.scale_yx
        return shapely.box(x0 * sx, y0 * sy, x1 * sx, y1 * sy)


def _to_tissue_bright_gray(img, channel_axis=0, p_low=1.0, p_high=99.0):
    """Convert a 2D/3D image to [0, 1] grayscale with tissue bright."""
    arr = np.asarray(img)
    if arr.ndim == 3:
        arr = arr.max(axis=channel_axis)
    elif arr.ndim != 2:
        raise ValueError(f"Expected a 2D or 3D image, got shape {arr.shape}")

    arr = arr.astype(np.float64)
    lo, hi = np.percentile(arr, (p_low, p_high))
    if hi <= lo:
        hi = lo + 1.0
    gray = np.clip((arr - lo) / (hi - lo), 0.0, 1.0)
    return 1.0 - gray if _should_invert(gray) else gray


def _should_invert(gray):
    if min(gray.shape) < 4:
        return False
    if float(np.percentile(gray, 95) - np.percentile(gray, 5)) < 0.15:
        return False

    bw = max(1, min(64, min(gray.shape) // 32))
    border = np.concatenate([
        gray[:bw].ravel(),
        gray[-bw:].ravel(),
        gray[:, :bw].ravel(),
        gray[:, -bw:].ravel(),
    ])
    return float(np.median(border)) > 0.65 and float(np.percentile(gray, 5)) < 0.55


def _tissue_mask(gray, policy=TissuePolicy()):
    """Build a cleaned tissue mask from an Otsu seed plus connected low-signal support."""
    threshold = skimage.filters.threshold_otsu(gray)
    seed = gray > threshold

    if policy.close_radius_px:
        seed = skimage.morphology.closing(
            seed, skimage.morphology.disk(policy.close_radius_px)
        )
    if policy.open_radius_px:
        seed = skimage.morphology.opening(
            seed, skimage.morphology.disk(policy.open_radius_px)
        )

    return _support_connected(gray, seed, threshold, policy)


def _support_connected(gray, seed, seed_threshold, policy):
    src = gray
    if policy.support_blur_sigma_px:
        src = skimage.filters.gaussian(
            gray,
            sigma=policy.support_blur_sigma_px,
            preserve_range=True,
        )

    bg_med, bg_mad = _border_background(gray)
    threshold = max(
        seed_threshold * policy.support_threshold_fraction,
        bg_med + policy.support_bg_sigmas * bg_mad,
        1e-6,
    )

    support = src > threshold
    if policy.support_density_threshold > 0:
        low_signal = gray > threshold
        density = skimage.filters.gaussian(
            low_signal.astype(float),
            sigma=policy.support_blur_sigma_px,
            preserve_range=True,
        )
        support |= density > policy.support_density_threshold

    support |= seed
    if policy.support_close_radius_px:
        support = skimage.morphology.closing(
            support, skimage.morphology.disk(policy.support_close_radius_px)
        )

    hole_size = max(64, int(np.ceil(policy.min_area_fraction * support.size)))
    support = skimage.morphology.remove_small_holes(support, max_size=hole_size)

    labels = skimage.measure.label(support)
    keep = np.unique(labels[seed])
    keep = keep[keep != 0]
    if len(keep) == 0:
        return np.zeros_like(seed, dtype=bool)
    return np.isin(labels, keep)


def _border_background(gray):
    bw = max(1, min(64, min(gray.shape) // 32))
    border = np.concatenate([
        gray[:bw].ravel(),
        gray[-bw:].ravel(),
        gray[:, :bw].ravel(),
        gray[:, -bw:].ravel(),
    ])
    median = float(np.median(border))
    mad = float(np.median(np.abs(border - median))) * 1.4826
    return median, mad


def _expand_bbox(bbox, shape, margin_fraction):
    y0, x0, y1, x1 = bbox
    my = int(np.ceil((y1 - y0) * margin_fraction))
    mx = int(np.ceil((x1 - x0) * margin_fraction))
    return (
        max(0, y0 - my),
        max(0, x0 - mx),
        min(shape[0], y1 + my),
        min(shape[1], x1 + mx),
    )


def _detect_tissue_boxes(gray, policy=TissuePolicy(), scale_yx=(1.0, 1.0)):
    """Detect tissue components and return largest-first ``TissueBox`` objects.

    With ``policy.attach_stray_radius_px`` > 0, each piece's box also grows to
    swallow small detached fragments (e.g. a severed thin tip) within that many
    smallest-level pixels. Only non-kept tissue ("strays") can be claimed, and
    each stray is assigned to its single nearest piece, so separate detected
    pieces are never merged.
    """
    mask = _tissue_mask(gray, policy)
    labels = skimage.measure.label(mask)
    min_area = policy.min_area_fraction * mask.size

    regions = [p for p in skimage.measure.regionprops(labels) if p.area >= min_area]
    regions.sort(key=lambda p: p.area, reverse=True)
    regions = regions[: policy.max_candidates]

    claim_label = _stray_claims(mask, labels, regions, policy.attach_stray_radius_px)

    boxes = []
    for region in regions:
        piece = labels == region.label
        if claim_label is not None:
            piece = piece | (claim_label == region.label)

        ys, xs = np.nonzero(piece)
        y0, x0, y1, x1 = _expand_bbox(
            (ys.min(), xs.min(), ys.max() + 1, xs.max() + 1),
            gray.shape,
            policy.box_margin_fraction,
        )
        component = piece[y0:y1, x0:x1]
        boxes.append(
            TissueBox(
                bbox_thumb=(y0, x0, y1, x1),
                mask=component,
                scale_yx=scale_yx,
                area_px=int(region.area),
            )
        )
    return boxes


def _stray_claims(mask, labels, regions, radius_px):
    """Assign each dropped fragment to its nearest kept piece within ``radius_px``.

    Returns an int label image where each stray pixel holds the label of the
    piece that claims it (0 elsewhere), or None when attachment is disabled.
    A single distance transform handles all pieces at once and guarantees a
    stray is claimed by only one (the nearest) piece.
    """
    radius = int(radius_px or 0)
    if radius <= 0 or not regions:
        return None

    kept_labels = [region.label for region in regions]
    kept_mask = np.isin(labels, kept_labels)
    strays = mask & ~kept_mask
    if not strays.any():
        return None

    distance, (iy, ix) = ndimage.distance_transform_edt(~kept_mask, return_indices=True)
    nearest = labels[iy, ix]  # label of the nearest kept component
    return np.where(strays & (distance <= radius), nearest, 0)


def detect_tissue(sdata, image_key="morphology_focus", policy=TissuePolicy()):
    """Detect on the smallest pyramid level; box geometry maps to level-0 pixels."""
    levels = pyramid_levels(sdata, image_key)
    smallest = levels[0]
    level0 = levels[-1]

    scale_yx = (
        int(level0.sizes["y"]) / int(smallest.sizes["y"]),
        int(level0.sizes["x"]) / int(smallest.sizes["x"]),
    )
    gray = _to_tissue_bright_gray(np.asarray(smallest))
    return _detect_tissue_boxes(gray, policy, scale_yx=scale_yx)


def _merge_group(boxes, order):
    """Union the 0-based indices in ``order`` into a single ``TissueBox``."""
    group = [boxes[i] for i in order]

    y0 = min(b.bbox_thumb[0] for b in group)
    x0 = min(b.bbox_thumb[1] for b in group)
    y1 = max(b.bbox_thumb[2] for b in group)
    x1 = max(b.bbox_thumb[3] for b in group)

    mask = np.zeros((y1 - y0, x1 - x0), dtype=bool)
    for box in group:
        by0, bx0, by1, bx1 = box.bbox_thumb
        mask[by0 - y0 : by1 - y0, bx0 - x0 : bx1 - x0] |= box.mask

    return TissueBox(
        bbox_thumb=(y0, x0, y1, x1),
        mask=mask,
        scale_yx=group[0].scale_yx,
        area_px=sum(box.area_px for box in group),
    )


def _normalize_merge_groups(indices):
    """Accept a flat list (one group) or a list of lists (several groups)."""
    seq = list(indices)
    if not seq:
        raise ValueError("merge_boxes needs at least one group of indices")

    nested = [isinstance(g, (list, tuple, set)) for g in seq]
    if all(nested):
        return [list(g) for g in seq]
    if not any(nested):
        return [seq]
    raise ValueError(
        "indices must be either a flat list of ints (one group) or a list of "
        "lists (several groups), not a mix of the two"
    )


def merge_boxes(boxes, indices):
    """Merge one or more groups of detected boxes (1-based, as numbered in the overview).

    ``indices`` may be either:

    * a flat list -- a single group, e.g. ``[2, 3]`` merges boxes 2 and 3; or
    * a list of lists -- several independent groups applied in one pass, e.g.
      ``[[2, 3], [4, 5]]`` merges 2+3 into one box and 4+5 into another.

    All groups are resolved against the original numbering simultaneously, so you
    don't need to account for renumbering between merges. Each group needs at
    least two distinct indices, and no index may appear in more than one group.

    For each group the chosen boxes' masks are unioned into their combined
    bounding box; boxes not named in any group are left unchanged. Use this to
    treat separate tissue pieces as a single crop when they are too far apart to
    merge via ``support_close_radius_px`` -- it is targeted (exactly the pieces
    you name) and distance-independent, so it never distorts outlines or pulls in
    unrelated pieces.

    Returns a new largest-first list; renumber by position in that list (so after
    merging, re-run ``preview_detection`` to see the new numbering).
    """
    groups = _normalize_merge_groups(indices)

    orders = []
    seen = {}
    for group in groups:
        order = sorted({i - 1 for i in group})
        if len(order) < 2:
            raise ValueError(
                f"each merge group needs at least two distinct indices; "
                f"group {group!r} does not"
            )
        if order[0] < 0 or order[-1] >= len(boxes):
            raise ValueError(f"indices {group} out of range for {len(boxes)} boxes")
        for k in order:
            if k in seen:
                raise ValueError(
                    f"box {k + 1} appears in more than one merge group "
                    f"({seen[k]!r} and {group!r})"
                )
            seen[k] = group
        orders.append(order)

    merged = [_merge_group(boxes, order) for order in orders]

    consumed = set(seen)
    remaining = [box for k, box in enumerate(boxes) if k not in consumed]
    result = remaining + merged
    result.sort(key=lambda box: box.area_px, reverse=True)
    return result
