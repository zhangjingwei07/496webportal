import math
from typing import Optional, Sequence

import plotly.graph_objects as go
from anndata import AnnData
from plotly.subplots import make_subplots

from .._utils import *


def rank_marker_genes(adata: AnnData,
                      key: str,
                      n_genes: int = 20,
                      gene_symbols: Optional[str] = None,
                      groups: Optional[Sequence[str]] = None,
                      ret_type: Optional[str] = 'json',
                      save: Optional[str] = None
                      ):
    data = adata.uns.get(key, None)
    if data is None:
        raise AttributeError("No groups found")
    reference = data['params']['reference']
    group_names = data['names'].dtype.names if groups is None else groups
    cols = math.ceil(len(group_names) ** 0.5)
    rows = math.ceil(len(group_names) / cols)

    fig = make_subplots(
        cols=cols,
        rows=rows,
        shared_yaxes=True,
        subplot_titles=[f" cluster {gn} vs {reference}" for gn in group_names])

    score_min, score_max = 0, 0
    for count, group_name in enumerate(group_names):
        gene_names = data['names'][group_name][:n_genes]
        scores = data['scores'][group_name][:n_genes]
        score_min = score_min if score_min <= scores.min() else scores.min()
        score_max = score_max if score_max >= scores.max() else scores.max()
        if gene_symbols is not None:
            if adata.raw is not None and data['params']['use_raw']:
                gene_names = [adata.raw.var[gene_symbols][g_name] for g_name in gene_names]
            else:
                gene_names = [adata.var[gene_symbols][g_name] for g_name in gene_names]
        fig.add_trace(go.Scatter(
            x=gene_names,
            y=scores,
            mode="markers",
            hovertemplate="gene: %{x}<br>score: %{y}",
            name=f" cluster {group_name}"
        ), row=count // cols + 1, col=count % cols + 1)
    fig.update_layout(showlegend=False, title="Rank Marker Genes")
    fig.update_yaxes(range=[score_min, score_max * 1.1])

    return fig_write_return(fig, ret_type=ret_type, save=save)
