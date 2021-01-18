from typing import Optional, Sequence

import numpy as np
import plotly.graph_objects as go
from anndata import AnnData
from .._utils import *

def tsne(adata: AnnData,
         names: Optional[Sequence[str]] = None,
         ret_type: Optional[str] = 'json',
         save: Optional[str] = None
         ):
    if names:
        clusterings = [n for n in names if n in adata.obs.columns]
        if clusterings:
            return _scatter_cluster(adata, "tSNE", clusterings, ret_type=ret_type, save=save)
    return _scatter(adata, "tSNE", names, ret_type=ret_type, save=save)


def umap(adata: AnnData,
         names: Optional[Sequence[str]] = None,
         ret_type: Optional[str] = 'json',
         save: Optional[str] = None
         ):
    if names:
        clusterings = [n for n in names if n in adata.obs.columns]
        if clusterings:
            return _scatter_cluster(adata, "UMAP", clusterings, ret_type=ret_type, save=save)
    return _scatter(adata, "UMAP", names, ret_type=ret_type, save=save)


def pca(adata: AnnData,
        names: Optional[Sequence[str]] = None,
        ret_type: Optional[str] = 'json',
        save: Optional[str] = None
        ):
    if names:
        clusterings = [n for n in names if n in adata.obs.columns]
        if clusterings:
            return _scatter_cluster(adata, "PCA", clusterings, ret_type=ret_type, save=save)
    return _scatter(adata, "PCA", names, ret_type=ret_type, save=save)


def _scatter_cluster(adata: AnnData,
                     basis: str,
                     clusterings: Sequence[str],
                     ret_type: Optional[str] = 'json',
                     save: Optional[str] = None
                     ) -> str:
    fig = go.Figure()
    visibility = []
    visible = True
    for method in clusterings:
        for cluster in adata.obs[method].cat.categories:
            fig.add_trace(
                go.Scattergl(
                    x=adata.obsm[f'X_{basis.lower()}'][adata.obs[method] == cluster, 0],
                    y=adata.obsm[f'X_{basis.lower()}'][adata.obs[method] == cluster, 1],
                    name=cluster,
                    mode='markers',
                    hoverinfo="text",
                    hovertext=adata.obs[adata.obs[method] == cluster].index,
                    text=np.where(adata.obs[method] == cluster)[0],
                    visible=visible
                )
            )
            visibility.append(method)
        visible = False
    buttons = []
    for method in clusterings:
        buttons.append({
            'args': [{'visible': [x == method for x in visibility]}, {'title': f'{basis} Plot - {method}'}],
            'label': method,
            'method': 'update'
        })
    _scatter_layout(fig, basis, buttons)

    return fig_write_return(fig, ret_type, save)


def _scatter(adata: AnnData,
             basis: str,
             names: Optional[Sequence[str]] = None,
             ret_type: Optional[str] = None,
             save: Optional[str] = None
             ) -> str:
    """
    Plot a 2D scatter plot on the given basis,
    :param adata: given annData that contains the obsm of the given basis
    :param basis: the name of the basis, currently accepting tsne and umap
    :param names: a list of gene names that will decide the coloring. If is None, chose the top 5 most.
    :param layer: Name of the AnnData object layer that wants to be plotted.
    :param use_raw: Use .raw attribute of adata for coloring with gene expression.
    """

    if f'X_{basis.lower()}' not in adata.obsm_keys():
        raise ValueError(f"Cannot find the basis. Was passed {basis}")

    if names is not None:
        names = names[:10]
    else:
        raise ValueError("No colors provided")

    matrix = adata[:, names].X
    if len(names) > 1:
        matrix = matrix.toarray()
    else:
        matrix = matrix[:, np.newaxis]

    fig = go.Figure()

    fig.add_trace(
        go.Scattergl(
            x=adata.obsm[f'X_{basis.lower()}'][:, 0],
            y=adata.obsm[f'X_{basis.lower()}'][:, 1],
            mode='markers',
            marker={
                'color': matrix[:, 0],
                'showscale': True,
                'opacity': .7,
                'colorbar': {
                    'outlinewidth': 0
                }
            },
            hovertemplate='%{hovertext}: %{marker.color}<extra></extra>',
            hovertext=adata.obs.index,
            text=np.arange(adata.n_obs)
        )
    )

    buttons = []
    for i in range(len(names)):
        buttons.append({
            'args': ['marker.color', [matrix[:, i]]],
            'label': names[i],
            'method': 'restyle'
        })

    _scatter_layout(fig, basis, buttons)

    return fig_write_return(fig, ret_type, save)


def _scatter_layout(fig, basis, buttons):
    fig.update_layout(
        template='plotly_dark',
        updatemenus=[
            go.layout.Updatemenu(buttons=buttons,
                                 showactive=True, active=0,
                                 direction="down", x=1, xanchor="right", y=1, yanchor="top"
                                 ),
        ],
        yaxis={'scaleanchor': "x", 'scaleratio': 1,
               'title': f'{basis}2',
               'showticklabels': False,
               'showgrid': False,
               'zeroline': False},
        xaxis={'title': f'{basis}1',
               'showticklabels': False,
               'showgrid': False,
               'zeroline': False},
        title=f"{basis} Plot - " + buttons[0]['label']
    )
