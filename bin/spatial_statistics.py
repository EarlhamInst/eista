#!/usr/bin/env python

# Set numba chache dir to current working directory (which is a writable mount also in containers)
import os
os.environ["NUMBA_CACHE_DIR"] = "."
os.environ[ 'MPLCONFIGDIR' ] = '/tmp/'

import scanpy as sc
import squidpy as sq
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import argparse
import sys, json
from pathlib import Path
from copy import deepcopy
import util

logger = util.get_named_logger('SPATIAL_STATISTICS')


def parse_args(argv=None):
    """Define and immediately parse command line arguments."""

    parser = argparse.ArgumentParser(
        description="Perform cell clustering and plot UMAPs of clustering",
    )
    parser.add_argument(
        "--h5ad",
        metavar="FILE_H5AD",
        type=Path,
        help="Input anndata data file.",
        required=True,
    )
    parser.add_argument(
        "--outdir",
        metavar="OUT_DIR",
        type=Path,
        help="Output directory.",
        required=True,
    )
    parser.add_argument(
        "--cluster_keys",
        default=None,
        help="Key in anndata.AnnData.obs where clustering is stored.",
    )
    parser.add_argument(
        "--coord_type",
        default='generic',
        choices=['grid', 'generic'],
        help="Type of coordinate system",
    )
    # parser.add_argument(
    #     "--autocorr_mode ",
    #     default='moran',
    #     choices=['moran', 'geary'],
    #     help="moran is for Moran’s I autocorrelation, and geary is for Geary’s C autocorrelation.",
    # )
    parser.add_argument(
        "--autocorr_n_perms",
        type=int,
        help="Number of permutations for the permutation test.",
        default=100,
    )
    parser.add_argument(
        "--autocorr_n_jobs",
        type=int,
        help="Number of parallel jobs.",
        default=5,
    )
    parser.add_argument(
        "--subsample_frac",
        type=float,
        help="Subsample to this fraction of the number of observations.",
        default=1,
    )
    # parser.add_argument(
    #     "--mode ",
    #     default='G',
    #     choices=['F', 'G', 'L'],
    #     help="Type of coordinate system",
    # )        
    # parser.add_argument(
    #     "--resolutions",
    #     type=util.floatlist,
    #     help="Resolution is used to control number of clusters.",
    #     default=None,
    # )
    parser.add_argument(
        "--meta",
        default='auto',
        choices=['auto', 'sample', 'group', 'plate'],
        help="Choose a metadata column as the batch for clustering",
    )
    parser.add_argument(
        "--fontsize",
        type=int,
        help="Set font size for plots.",
        default=12,
    )
    parser.add_argument(
        "--pdf",
        help="Whether to generate figure files in PDF format.",
        action='store_true',
    )                      
    return parser.parse_args(argv)


def main(argv=None):

    args = parse_args(argv)

    if not args.h5ad.is_file():
        logger.error(f"The given input file {args.h5ad} was not found!")
        sys.exit(2)

    plt.rcParams.update({
        "font.size": args.fontsize,
        # "axes.titlesize": 'medium',
        # "axes.labelsize": 'small',
        # "xtick.labelsize": 'small',
        # "ytick.labelsize": 'small',
        # "legend.fontsize": 'small',
    })

    util.check_and_create_folder(args.outdir)
    path_statistics = Path(args.outdir)
    util.check_and_create_folder(path_statistics)

    adata = sc.read_h5ad(args.h5ad)

    if args.meta == 'auto':
        # batch = 'group' if hasattr(adata.obs, 'group') else 'sample'
        batch = 'sample'
        if hasattr(adata.obs, 'group'):
            batch = 'group'
        elif hasattr(adata.obs, 'plate'):
            batch = 'plate' 
    else:
        batch = args.meta
    
    
    # Building the spatial neighbors graphs
    # sq.gr.spatial_neighbors(adata, coord_type=args.coord_type, delaunay=True)


    cluster_keys = []
    if not args.cluster_keys:
        cluster_keys = [col for col in adata.obs.columns if col.startswith('leiden')]
    else:
        cluster_keys = args.cluster_keys.split(',')

    for sid in sorted(adata.obs[batch].unique()):
        adata_s = adata[adata.obs[batch]==sid]   
        path_statistics_s = Path(path_statistics, f"{batch}_{sid}")
        util.check_and_create_folder(path_statistics_s)
        sq.gr.spatial_neighbors(adata_s, coord_type=args.coord_type, delaunay=True)
        for cluster_key in cluster_keys:
            ser_counts = adata_s.obs[cluster_key].value_counts()
            ser_counts.name = "cell counts"
            meta_leiden = pd.DataFrame(ser_counts)

            # Compute centrality scores
            sq.gr.centrality_scores(adata_s, cluster_key=cluster_key)
            df_central = deepcopy(adata_s.uns[f"{cluster_key}_centrality_scores"])
            df_central.index = meta_leiden.index.tolist()
            ser_closeness = df_central["closeness_centrality"].sort_values(ascending=False)
            ser_degree = df_central["degree_centrality"].sort_values(ascending=False)
            ser_cluster = df_central["average_clustering"].sort_values(ascending=False)            
            with plt.rc_context():
                sq.pl.centrality_scores(
                    adata_s, 
                    cluster_key=cluster_key, 
                    figsize=(16, 5)
                )
                plt.savefig(Path(path_statistics_s, f"centrality_scores_{cluster_key}.png"), bbox_inches="tight")
                if args.pdf:
                    plt.savefig(Path(path_statistics_s, f"centrality_scores_{cluster_key}.pdf"), bbox_inches="tight")
            with plt.rc_context():
                sq.pl.spatial_scatter(
                    adata_s, 
                    groups=ser_closeness.index.tolist()[:min(5, ser_closeness.size//2)], 
                    shape=None, 
                    color=cluster_key,
                    wspace=0.4,
                    title=f"High closeness - {cluster_key}",
                )
                plt.savefig(Path(path_statistics_s, f"spatial_scatter_closeness_high_{cluster_key}.png"), bbox_inches="tight")
                if args.pdf:
                    plt.savefig(Path(path_statistics_s, f"spatial_scatter_closeness_high_{cluster_key}.pdf"), bbox_inches="tight")
            with plt.rc_context():
                sq.pl.spatial_scatter(
                    adata_s, 
                    groups=ser_closeness.index.tolist()[max(-5, -ser_closeness.size//2):], 
                    shape=None, 
                    color=cluster_key,
                    wspace=0.4,
                    title=f"Low closeness - {cluster_key}",
                )
                plt.savefig(Path(path_statistics_s, f"spatial_scatter_closeness_low_{cluster_key}.png"), bbox_inches="tight")
                if args.pdf:
                    plt.savefig(Path(path_statistics_s, f"spatial_scatter_closeness_low_{cluster_key}.pdf"), bbox_inches="tight")
            with plt.rc_context():
                sq.pl.spatial_scatter(
                    adata_s, 
                    groups=ser_degree.index.tolist()[:min(5, ser_degree.size//2)], 
                    shape=None, 
                    color=cluster_key,
                    wspace=0.4,
                    title=f"High degree - {cluster_key}",
                )
                plt.savefig(Path(path_statistics_s, f"spatial_scatter_degree_high_{cluster_key}.png"), bbox_inches="tight")
                if args.pdf:
                    plt.savefig(Path(path_statistics_s, f"spatial_scatter_degree_high_{cluster_key}.pdf"), bbox_inches="tight")
            with plt.rc_context():
                sq.pl.spatial_scatter(
                    adata_s, 
                    groups=ser_degree.index.tolist()[max(-5, -ser_degree.size//2):], 
                    shape=None, 
                    color=cluster_key,
                    wspace=0.4,
                    title=f"Low degree - {cluster_key}",
                )
                plt.savefig(Path(path_statistics_s, f"spatial_scatter_degree_low_{cluster_key}.png"), bbox_inches="tight")
                if args.pdf:
                    plt.savefig(Path(path_statistics_s, f"spatial_scatter_degree_low_{cluster_key}.pdf"), bbox_inches="tight")
            with plt.rc_context():
                sq.pl.spatial_scatter(
                    adata_s, 
                    groups=ser_cluster.index.tolist()[:min(5, ser_cluster.size//2)], 
                    shape=None, 
                    color=cluster_key,
                    wspace=0.4,
                    title=f"High clustering coefficient - {cluster_key}",
                )
                plt.savefig(Path(path_statistics_s, f"spatial_scatter_clustering_high_{cluster_key}.png"), bbox_inches="tight")
                if args.pdf:
                    plt.savefig(Path(path_statistics_s, f"spatial_scatter_clustering_high_{cluster_key}.pdf"), bbox_inches="tight")
            with plt.rc_context():
                sq.pl.spatial_scatter(
                    adata_s, 
                    groups=ser_cluster.index.tolist()[max(-5, -ser_cluster.size//2):], 
                    shape=None, 
                    color=cluster_key,
                    wspace=0.4,
                    title=f"Low clustering coefficient - {cluster_key}",
                )
                plt.savefig(Path(path_statistics_s, f"spatial_scatter_clustering_low_{cluster_key}.png"), bbox_inches="tight")
                if args.pdf:
                    plt.savefig(Path(path_statistics_s, f"spatial_scatter_clustering_low_{cluster_key}.pdf"), bbox_inches="tight")

            # Neighbors enrichment analysis
            sq.gr.nhood_enrichment(adata_s, cluster_key=cluster_key, show_progress_bar=False)
            with plt.rc_context():
                sq.pl.nhood_enrichment(
                    adata_s,
                    cluster_key=cluster_key,
                    figsize=(8, 8),
                    title="Neighborhood enrichment heatmap",
                )
                plt.savefig(Path(path_statistics_s, f"neighbors_enrichment_{cluster_key}.png"), bbox_inches="tight")
                if args.pdf:
                    plt.savefig(Path(path_statistics_s, f"neighbors_enrichment_{cluster_key}.pdf"), bbox_inches="tight")
    
            # Ripley’s statistics
            # sq.gr.ripley(adata_s, cluster_key=cluster_key, mode='L')
            # with plt.rc_context():
            #     sq.pl.ripley(
            #         adata_s, 
            #         cluster_key="leiden", 
            #         mode='L',
            #     )
            #     plt.savefig(Path(path_statistics_s, f"Ripley_L_{cluster_key}.png"), bbox_inches="tight")
            # sq.gr.ripley(adata_s, cluster_key=cluster_key, mode='F')
            # with plt.rc_context():
            #     sq.pl.ripley(
            #         adata_s, 
            #         cluster_key="leiden", 
            #         mode='F',
            #     )
            #     plt.savefig(Path(path_statistics_s, f"Ripley_F_{cluster_key}.png"), bbox_inches="tight")

        # Moran’s I score
        adata_subsample = adata_s
        if args.subsample_frac < 1:
            adata_subsample = sc.pp.subsample(adata_s, fraction=args.subsample_frac, copy=True)
            sq.gr.spatial_neighbors(adata_subsample, coord_type=args.coord_type, delaunay=True)
        sq.gr.spatial_autocorr(
            adata_subsample,
            mode='moran',
            n_perms=args.autocorr_n_perms,
            n_jobs=args.autocorr_n_jobs,
        )
        # amode = {'moran': 'moranI', 'geary': 'gearyC'}.get(args.autocorr_mode)
        pd.concat([adata_subsample.uns['moranI'].head(50), adata_subsample.uns['moranI'].tail(50)]) \
            .reset_index().rename(columns={'index': 'Gene'}) \
            .to_csv(Path(path_statistics_s, 'autocorr_moranI.csv'), index=False)
        # adata_subsample.uns['moran'].head(100).to_csv(Path(path_statistics_s, 'autocorr_moranI.csv'), index=False)
        for gene in adata_subsample.uns['moranI'].head(6).index.tolist():
            with plt.rc_context():
                sq.pl.spatial_scatter(
                    adata_subsample,
                    color=[gene],
                    shape=None,
                    size=2,
                    img=False,
                )
                plt.savefig(Path(path_statistics_s, f"spatial_scatter_top_{gene}.png"), bbox_inches="tight")
                if args.pdf:
                    plt.savefig(Path(path_statistics_s, f"spatial_scatter_top_{gene}.pdf"), bbox_inches="tight")
        for gene in adata_subsample.uns['moranI'].tail(6).index.tolist():
            with plt.rc_context():
                sq.pl.spatial_scatter(
                    adata_subsample,
                    color=[gene],
                    shape=None,
                    size=2,
                    img=False,
                )
                plt.savefig(Path(path_statistics_s, f"spatial_scatter_bot_{gene}.png"), bbox_inches="tight")
                if args.pdf:
                    plt.savefig(Path(path_statistics_s, f"spatial_scatter_bot_{gene}.pdf"), bbox_inches="tight")



    # save analysis parameters into a json file
    with open(Path(path_statistics, 'parameters.json'), 'w') as file:
        params = {}
        params.update({"--h5ad": str(args.h5ad)})        
        params.update({"--coord_type": str(args.coord_type)})        
        if args.cluster_keys: params.update({"--cluster_keys": args.cluster_keys})        
        # if args.integrate: params.update({"--integrate": args.integrate})        
        params.update({"--autocorr_n_perms": args.autocorr_n_perms})        
        params.update({"--autocorr_n_jobs": args.autocorr_n_jobs})        
        params.update({"--subsample_frac": args.subsample_frac})        
        params.update({"--meta": args.meta})        
        json.dump(params, file, indent=4)



if __name__ == "__main__":
    sys.exit(main())
