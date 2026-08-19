# SpatialData Dissect ✂️

Tissue detection and cropping for [SpatialData](https://spatialdata.scverse.org/)
/ [Xenium](https://www.10xgenomics.com/platforms/xenium) images. `spatialdata_dissect` finds the main tissue pieces on the smallest pyramid
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

#### Parameters

`TissuePolicy` is a frozen dataclass; pass any subset of fields to override the
defaults. The knobs fall into four groups that follow the detection pipeline:

- **seed cleanup** — `close_radius_px`, `open_radius_px`
- **support expansion** (how much faint tissue is pulled in around the bright
  seed) — `support_blur_sigma_px`, `support_threshold_fraction`,
  `support_bg_sigmas`, `support_density_threshold`, `support_close_radius_px`
- **component selection** — `min_area_fraction`, `max_candidates`
- **box building** — `box_margin_fraction`, `attach_stray_radius_px`

The support-expansion group are the recall levers: they decide how much
low-density tissue (e.g. dermis) survives. The selection and box knobs only
affect how captured tissue is split into pieces and cropped, so they won't help
if tissue was never captured in the first place — tune the support group first.

| Parameter | Default | What it does |
| --- | --- | --- |
| `min_area_fraction` | `0.002` | Minimum component area, as a fraction of the image, for a blob to count as a real tissue piece; the same value sets the hole-filling size in the support step. Raise it to ignore small fragments and fill larger interior gaps; lower it to keep tiny pieces and preserve small holes. |
| `close_radius_px` | `4` | Radius of the morphological closing (dilate then erode) applied to the Otsu seed. Fills pin-holes and bridges nearly-touching bright pixels, so a speckled dense region becomes one solid blob. Raise it to consolidate a broken-up seed; lower it to keep fine structure separate. |
| `open_radius_px` | `2` | Radius of the opening (erode then dilate) applied right after closing; deletes isolated specks smaller than the disk. Raise it if the seed picks up scattered noise; set it to `0` to keep every bright speck. |
| `support_blur_sigma_px` | `4.0` | Gaussian blur applied to the grayscale before the support threshold (also the sigma used by the density path). Larger values recover sparser tissue but smear intensities down; `0` thresholds the raw image. |
| `support_threshold_fraction` | `0.15` | Sets the support threshold to at least this fraction of the Otsu seed threshold — the primary "how faint can included tissue be" knob. Lower includes fainter tissue; too low starts pulling in background. |
| `support_bg_sigmas` | `3.0` | Floors the support threshold at background median + this many robust standard deviations, keeping it from sinking into background noise. Usually inert on clean black backgrounds; on noisy ones it prevents flooding. Raise it to be more conservative near noise; lower it to allow the threshold closer to background. |
| `support_density_threshold` | `0.05` | When positive, counts the local fraction of above-threshold pixels and includes any neighbourhood where that fraction exceeds this value, even if no single pixel is bright (a fraction between 0 and 1). This is the lever for sparse-but-real tissue such as dermis; `0` disables it. |
| `support_close_radius_px` | `8` | Radius of the closing applied to the combined support mask, joining pieces that are close. Larger bridges more distant fragments; too large can fuse things that should stay separate. |
| `box_margin_fraction` | `0.05` | Pads each piece's bounding box by this fraction of its own width and height on every side. Affects only the crop rectangle, not the tissue mask. Raise it if crops are clipping tissue at the edges. |
| `attach_stray_radius_px` | `25` | Grows each piece's box to swallow detached fragments within this many smallest-level pixels; `0` disables. Only claims strays (non-kept tissue), so it never merges two separate detected pieces. |
| `max_candidates` | `12` | Caps how many pieces are returned, after sorting largest-first. With more physical sections than this you'll lose the smallest; lower it to trim spurious extras to the biggest few. |

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
