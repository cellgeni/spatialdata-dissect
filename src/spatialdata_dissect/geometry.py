"""Convert detected boxes into Shapely geometry and map to a coordinate system.

``get_dissected_boxes`` is the one-call convenience wrapper: it turns a list of
detected boxes straight into query-ready polygons, choosing how overlapping
neighbours are handled and which coordinate system the result lives in.
"""

import numpy as np
import shapely
import skimage.measure
import spatialdata


def _other_tissue_mask(box, boxes):
    """Other detected components, clipped into ``box`` at thumbnail resolution."""
    y0, x0, y1, x1 = box.bbox_thumb
    out = np.zeros((y1 - y0, x1 - x0), dtype=bool)

    for other in boxes:
        if other is box:
            continue

        oy0, ox0, oy1, ox1 = other.bbox_thumb
        iy0, ix0 = max(y0, oy0), max(x0, ox0)
        iy1, ix1 = min(y1, oy1), min(x1, ox1)
        if iy0 >= iy1 or ix0 >= ix1:
            continue

        out[iy0 - y0 : iy1 - y0, ix0 - x0 : ix1 - x0] |= other.mask[
            iy0 - oy0 : iy1 - oy0,
            ix0 - ox0 : ix1 - ox0,
        ]

    return out


def _mask_geometry(mask, offset_yx, scale_yx, simplify_tol=1.5):
    """Convert disconnected mask components to a Shapely geometry."""
    labels = skimage.measure.label(mask)
    polygons = []
    y0, x0 = offset_yx
    sy, sx = scale_yx

    for label_id in range(1, labels.max() + 1):
        component = labels == label_id
        contours = skimage.measure.find_contours(
            np.pad(component, 1).astype(float), 0.5
        )
        if not contours:
            continue

        contour = max(contours, key=len) - 1
        xy = np.column_stack([
            (contour[:, 1] + x0) * sx,
            (contour[:, 0] + y0) * sy,
        ])
        polygon = shapely.Polygon(xy)
        if not polygon.is_valid:
            polygon = polygon.buffer(0)
        if simplify_tol and not polygon.is_empty:
            polygon = polygon.simplify(simplify_tol)
        if not polygon.is_empty:
            polygons.append(polygon)

    if not polygons:
        return None
    return shapely.union_all(polygons)


def rect_minus_others(
    box,
    boxes,
    pad_others=0.0,
    simplify_tol=1.5,
    single=True,
):
    """Return this full bounding box minus tissue belonging to the other boxes.

    The result is in level-0 pixel coordinates and uses the same subtraction mask as ``preview_detection``.
    """
    rect = box.rect_polygon()
    y0, x0, _, _ = box.bbox_thumb
    others = _other_tissue_mask(box, boxes)

    geometry = _mask_geometry(
        others,
        offset_yx=(y0, x0),
        scale_yx=box.scale_yx,
        simplify_tol=simplify_tol,
    )
    if geometry is not None:
        if pad_others:
            geometry = geometry.buffer(pad_others)
        rect = rect.difference(geometry)

    if simplify_tol and not rect.is_empty:
        rect = rect.simplify(simplify_tol)
    if single and rect.geom_type == "MultiPolygon" and len(rect.geoms):
        rect = max(rect.geoms, key=lambda geom: geom.area)
    return rect


def to_coordinate_system(sdata, poly, coordinate_system="global", image_key="morphology_focus"):
    """Map level-0 pixel geometry into a named SpatialData coordinate system."""
    transform = spatialdata.transformations.get_transformation(
        sdata.images[image_key],
        to_coordinate_system=coordinate_system,
    )
    matrix = transform.to_affine_matrix(
        input_axes=("x", "y"),
        output_axes=("x", "y"),
    )

    def apply(points):
        homogeneous = np.column_stack([points, np.ones(len(points))])
        return (matrix @ homogeneous.T).T[:, :2]

    return shapely.transform(poly, apply)


def to_global_coordinates(sdata, poly, image_key="morphology_focus"):
    """Map level-0 pixel geometry to the SpatialData ``global`` coordinate system.

    Thin wrapper around ``to_coordinate_system(..., coordinate_system="global")``.
    """
    return to_coordinate_system(sdata, poly, coordinate_system="global", image_key=image_key)


_KEEP_ALIASES = {None, "keep", "none", "full", "rect", "rectangle"}
_MINUS_ALIASES = {"minus", "subtract", "difference"}
_PIXEL_ALIASES = {None, "pixel", "pixels", "none"}


def get_dissected_boxes(
    sdata,
    boxes,
    overlap="minus",
    transformation="global",
    image_key="morphology_focus",
    pad_others=0.0,
    simplify_tol=1.5,
    single=True,
):
    """Turn detected boxes into query-ready polygons in one call.

    Equivalent to, but tidier than::

        [to_global_coordinates(sdata, rect_minus_others(b, boxes)) for b in boxes]

    Parameters
    ----------
    sdata, boxes:
        The SpatialData object and the ``TissueBox`` list from ``detect_tissue``.
    overlap:
        How to treat tissue that belongs to *other* boxes but falls inside this
        box's rectangle. ``"minus"`` (default) subtracts it (see
        ``rect_minus_others``); ``"keep"`` (aliases: ``None``, ``"full"``,
        ``"rect"``) returns the plain bounding rectangle instead.
    transformation:
        Target coordinate system name, e.g. ``"global"`` (default), so the
        polygons can be passed straight to ``spatialdata.polygon_query``. Use
        ``None`` or ``"pixel"`` to leave them in level-0 pixel coordinates.
    pad_others, simplify_tol, single:
        Forwarded to ``rect_minus_others`` when ``overlap="minus"``.

    Returns
    -------
    list of shapely geometries, one per box, in the same order as ``boxes``.
    """
    key = overlap.lower() if isinstance(overlap, str) else overlap
    if key in _MINUS_ALIASES:
        make = lambda box: rect_minus_others(
            box, boxes, pad_others=pad_others, simplify_tol=simplify_tol, single=single
        )
    elif key in _KEEP_ALIASES:
        make = lambda box: box.rect_polygon()
    else:
        raise ValueError(f"overlap must be 'minus' or 'keep', got {overlap!r}")

    tkey = transformation.lower() if isinstance(transformation, str) else transformation
    in_pixels = tkey in _PIXEL_ALIASES

    polys = []
    for box in boxes:
        poly = make(box)
        if not in_pixels:
            poly = to_coordinate_system(
                sdata, poly, coordinate_system=transformation, image_key=image_key
            )
        polys.append(poly)
    return polys
