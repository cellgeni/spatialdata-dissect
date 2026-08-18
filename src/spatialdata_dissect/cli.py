"""Command-line interface for batch tissue detection and preview.

Two ways to point it at data:

* ``--csv experiments.csv`` -- a table with one row per experiment. By default
  it reads the ``region_name`` and ``xenium_exp`` columns (override with
  ``--region-column`` / ``--path-column``). Each region gets its own
  subdirectory under ``--outdir``.

* one or more positional Xenium experiment directories -- the directory's base
  name is used as the region name.

Examples::

    spatialdata_dissect --csv experiments.csv --outdir results/
    spatialdata_dissect /path/to/xenium_exp_a /path/to/xenium_exp_b --outdir results/
"""

import argparse
import os
import sys
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from .detection import TissuePolicy, detect_tissue
from .plot import preview_detection


def _build_parser():
    parser = argparse.ArgumentParser(
        prog="spatialdata_dissect",
        description="Detect tissue pieces in Xenium/SpatialData images and write detection previews.",
    )
    parser.add_argument(
        "xenium",
        nargs="*",
        help="Xenium experiment directories to process (ignored if --csv is given).",
    )
    parser.add_argument(
        "--csv",
        help="CSV listing experiments to process, one per row.",
    )
    parser.add_argument(
        "--region-column",
        default="region_name",
        help="CSV column holding the region name (default: region_name).",
    )
    parser.add_argument(
        "--path-column",
        default="xenium_exp",
        help="CSV column holding the Xenium experiment path (default: xenium_exp).",
    )
    parser.add_argument(
        "-o",
        "--outdir",
        default="tissue_detection",
        help="Output directory; each region is written to a subdirectory of it.",
    )
    parser.add_argument(
        "--image-key",
        default="morphology_focus",
        help="Image key to detect on (default: morphology_focus).",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=150,
        help="DPI for the rendered previews (default: 150).",
    )
    parser.add_argument(
        "--no-overview",
        dest="save_overview",
        action="store_false",
        help="Do not save the whole-image overview with bounding boxes.",
    )
    parser.add_argument(
        "--no-crops",
        dest="save_crops",
        action="store_false",
        help="Do not save the individual per-box crop previews.",
    )
    parser.add_argument(
        "--max-candidates",
        type=int,
        default=TissuePolicy.max_candidates,
        help="Maximum number of tissue pieces to keep per image.",
    )
    parser.add_argument(
        "--attach-stray-radius-px",
        type=int,
        default=TissuePolicy.attach_stray_radius_px,
        help="Grow each box to swallow detached fragments within this many smallest-level pixels (0 disables).",
    )
    return parser


def _iter_experiments(args):
    """Yield ``(region_name, xenium_path)`` pairs from the parsed arguments."""
    if args.csv:
        import pandas as pd

        df = pd.read_csv(args.csv)
        for _, row in df.iterrows():
            yield str(row[args.region_column]), str(row[args.path_column])
    else:
        for path in args.xenium:
            region = os.path.basename(os.path.normpath(path))
            yield region, path


def main(argv=None):
    warnings.filterwarnings("ignore", category=UserWarning)

    args = _build_parser().parse_args(argv)

    if not args.csv and not args.xenium:
        print(
            "error: provide --csv or at least one Xenium experiment directory.",
            file=sys.stderr,
        )
        return 2

    import spatialdata_io

    policy = TissuePolicy(
        max_candidates=args.max_candidates,
        attach_stray_radius_px=args.attach_stray_radius_px,
    )

    experiments = list(_iter_experiments(args))
    for region_name, xenium_path in experiments:
        print(region_name)
        sdata = spatialdata_io.xenium(xenium_path)
        boxes = detect_tissue(sdata, image_key=args.image_key, policy=policy)
        overview_fig, crop_figs = preview_detection(
            sdata,
            boxes,
            outdir=os.path.join(args.outdir, region_name),
            image_key=args.image_key,
            dpi=args.dpi,
            save_overview=args.save_overview,
            save_crops=args.save_crops,
        )
        for fig in filter(None, [overview_fig, *crop_figs]):
            plt.close(fig)

    print(f"Done: processed {len(experiments)} experiment(s).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
