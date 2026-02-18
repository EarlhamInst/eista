#!/usr/bin/env python

# Set numba chache dir to current working directory (which is a writable mount also in containers)
import os
os.environ["NUMBA_CACHE_DIR"] = "."
os.environ[ 'MPLCONFIGDIR' ] = '/tmp/'


import scanpy as sc
import squidpy as sq
import pandas as pd
import argparse
import sys
from pathlib import Path
import spatialdata as sd
import util

logger = util.get_named_logger('SPATIAL_TO_H5AD')


def parse_args(argv=None):
    """Define and immediately parse command line arguments."""

    parser = argparse.ArgumentParser(
        description="Read spatial data files and convert to h5ad file",
    )
    parser.add_argument(
        "--tech",
        help="The spatial single-cell technology.",
        required=True,
    )
    parser.add_argument(
        "--outfile",
        metavar="OUT_file",
        type=Path,
        help="Output file.",
        required=True,
    )
    parser.add_argument(
        "--datadir",
        metavar="DATA_DIR",
        type=Path,
        help="Raw data directory.",
        required=True,
    )
    parser.add_argument(
        "--counts",
        metavar="COUNTS_FILE",
        type=Path,
        help="The cell-by-gene counts file.",
    )
    parser.add_argument(
        "--metadata",
        metavar="META_FILE",
        type=Path,
        help="The cell metadata file.",
    )
    parser.add_argument(
        "--transformation",
        metavar="TRANSFORMATION_FILE",
        type=Path,
        help="The micron_to_mosaic_pixel_transform file.",
    )          
    return parser.parse_args(argv)


def main(argv=None):

    args = parse_args(argv)
    adata = None

    if args.tech == 'vizgen':
        if not args.counts.is_file():
            logger.error(f"The given input file {args.counts} was not found!")
            sys.exit(2)

        if args.metadata and not args.metadata.is_file():
            logger.error(f"The given input file {args.metadata} was not found!")
            sys.exit(2)

        # transformation = Path(args.datadir, "images/micron_to_mosaic_pixel_transform.csv")


        adata = sq.read.vizgen(
            path=args.datadir,
            counts_file=str(args.counts),
            meta_file=str(args.metadata),
            transformation_file=str(args.transformation) if args.transformation.is_file() else None,
        )
    elif args.tech == 'xenium':
        from spatialdata_io import xenium
        xenium_path = args.datadir
        sdata = xenium(xenium_path)
        adata = sdata.tables["table"]


    # save the AnnData into a h5ad file
    if adata:
        adata.write_h5ad(args.outfile, compression="gzip")
        # adata.write_h5ad(Path(args.outdir, 'spatial_converted.h5ad'))



if __name__ == "__main__":
    sys.exit(main())
