# SpatialData Dissect ✂️

Tissue detection and cropping for [SpatialData](https://spatialdata.scverse.org/)
/ Xenium images. `spatialdata_dissect` finds the main tissue pieces on the smallest pyramid
level of a morphology image, draws diagnostic previews, and produces
level-0-pixel or global-coordinate geometry you can feed straight into
`spatialdata.polygon_query`.

## Install

From the project repository:

```bash
pip install git+https://github.com/cellgeni/spatialdata-dissect
```

## Library usage

```python
import spatialdata
import spatialdata_io
import spatialdata_dissect

sdata = spatialdata_io.xenium("path/to/xenium_experiment")

# Detect tissue pieces (largest-first).
boxes = spatialdata_dissect.detect_tissue(sdata)

# In a notebook: build the figures in memory and display them inline. Nothing
# is written to disk unless you pass `outdir`.
overview_fig, crop_figs = spatialdata_dissect.preview_detection(sdata, boxes)

# Write an overview image plus one crop preview per detected piece to disk.
spatialdata_dissect.preview_detection(sdata, boxes, outdir="out")

# You can also build a single figure directly:
#   fig = spatialdata_dissect.overview_figure(np.asarray(smallest_level), boxes)
#   fig = spatialdata_dissect.crop_figure(sdata, boxes[0], boxes)

# Optionally merge separate pieces into a single crop. Indices are 1-based,
# matching the numbers drawn in the overview image. Re-run preview_detection
# afterwards to see the new numbering.
# boxes = spatialdata_dissect.merge_boxes(boxes, [2, 3])            # one group: merge 2+3
# boxes = spatialdata_dissect.merge_boxes(boxes, [[2, 3], [4, 5]])  # two groups in one pass

# Turn every box into a query-ready polygon in one call. `overlap` chooses
# whether neighbouring tissue is subtracted ("minus") or the plain rectangle is
# kept ("keep"); `transformation` is the target coordinate system (or None for
# level-0 pixels).
polys = spatialdata_dissect.get_dissected_boxes(sdata, boxes, overlap="minus", transformation="global")

for poly in polys:
    cropped = spatialdata.polygon_query(sdata, poly, "global", filter_table=True)
    # ... save or process `cropped`

# The lower-level pieces are still available if you need them:
#   poly = spatialdata_dissect.rect_minus_others(box, boxes)
#   poly = spatialdata_dissect.to_coordinate_system(sdata, poly, "global")
```

### Tuning detection

`detect_tissue` accepts a `TissuePolicy` with all detection knobs:

```python
from spatialdata_dissect import TissuePolicy, detect_tissue

policy = TissuePolicy(
    max_candidates=6,
    attach_stray_radius_px=25,   # swallow detached fragments; 0 disables
    support_close_radius_px=8,   # bridge nearby fragments into one piece
)
boxes = detect_tissue(sdata, policy=policy)
```

## Command-line usage

The package installs a `spatialdata_dissect` console script.

Process a batch described by a CSV with `region_name` and `xenium_exp` columns:

```bash
spatialdata_dissect --csv experiments.csv --outdir results/
```

Or point it directly at one or more Xenium experiment directories (the
directory base name becomes the region name):

```bash
spatialdata_dissect /path/to/xenium_exp_a /path/to/xenium_exp_b --outdir results/
```

Save the overview and the individual crops independently with `--no-overview`
and `--no-crops` (both are saved by default):

```bash
spatialdata_dissect --csv experiments.csv --outdir results/ --no-crops      # overview only
spatialdata_dissect --csv experiments.csv --outdir results/ --no-overview   # crops only
```

Other useful options: `--image-key`, `--dpi`, `--max-candidates`,
`--attach-stray-radius-px`, and `--region-column` / `--path-column` to override
the CSV column names. Run `spatialdata_dissect --help` for the full list.

## Public API

| Name | Purpose |
| --- | --- |
| `detect_tissue(sdata, ...)` | Detect tissue pieces, return largest-first `TissueBox` list. |
| `merge_boxes(boxes, indices)` | Fuse one group `[2,3]` or several `[[2,3],[4,5]]` (1-based) in one pass. |
| `get_dissected_boxes(sdata, boxes, overlap=, transformation=)` | One call: boxes → query-ready polygons. |
| `rect_minus_others(box, boxes, ...)` | Box rectangle minus other pieces' tissue (level-0 pixels). |
| `to_coordinate_system(sdata, poly, coordinate_system=, ...)` | Map pixel geometry to any named coordinate system. |
| `to_global_coordinates(sdata, poly, ...)` | Shorthand for `to_coordinate_system(..., "global")`. |
| `preview_detection(sdata, boxes, outdir=None, ...)` | Build `(overview_fig, crop_figs)`; save PNGs only if `outdir` is given. |
| `overview_figure(image, boxes, ...)` | Build just the overview figure. |
| `crop_figure(sdata, box, boxes, ...)` | Build just one crop figure. |
| `TissuePolicy` | Detection parameters. |
| `TissueBox` | A detected component and its bounding box. |

`preview_detection` returns the matplotlib figures it builds, so it works inline
in a notebook without touching disk; pass `outdir="..."` to also save PNGs. When
looping over many images, close the returned figures with `plt.close(fig)` to
free memory (the CLI does this for you).

## Package layout

```
src/spatialdata_dissect/
├── __init__.py      public API re-exports
├── detection.py     TissuePolicy, TissueBox, detect_tissue, merge_boxes + helpers
├── pyramid.py       pyramid level access
├── geometry.py      rect_minus_others, to_global_coordinates + helpers
├── plot.py          preview_detection, overview_figure, crop_figure
└── cli.py           command-line entry point
```
