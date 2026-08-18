"""spatialdata_dissect -- tissue detection and cropping for SpatialData/Xenium images.

Typical use::

    import spatialdata_io
    import spatialdata_dissect

    sdata = spatialdata_io.xenium("path/to/xenium_experiment")
    boxes = spatialdata_dissect.detect_tissue(sdata)

    # In a notebook: build figures in memory and display them inline, no files.
    overview_fig, crop_figs = spatialdata_dissect.preview_detection(sdata, boxes)

    # Or write PNGs to disk by passing an output directory.
    spatialdata_dissect.preview_detection(sdata, boxes, outdir="out")

    # Optionally treat separate pieces as one crop (1-based, as numbered in
    # the overview image), then re-preview to see the new numbering:
    # boxes = spatialdata_dissect.merge_boxes(boxes, [2, 3])

    # Turn a box into a global-coordinate polygon for cropping:
    for box in boxes:
        poly = spatialdata_dissect.rect_minus_others(box, boxes)
        query = spatialdata_dissect.to_global_coordinates(sdata, poly)
        # cropped = spatialdata.polygon_query(sdata, query, "global", filter_table=True)
"""

from .detection import TissueBox, TissuePolicy, detect_tissue, merge_boxes
from .geometry import (
    get_dissected_boxes,
    rect_minus_others,
    to_coordinate_system,
    to_global_coordinates,
)
from .plot import crop_figure, overview_figure, preview_detection

__all__ = [
    "TissuePolicy",
    "TissueBox",
    "detect_tissue",
    "merge_boxes",
    "get_dissected_boxes",
    "rect_minus_others",
    "to_coordinate_system",
    "to_global_coordinates",
    "preview_detection",
    "overview_figure",
    "crop_figure",
]

__version__ = "0.1.0"
