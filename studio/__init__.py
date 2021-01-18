from typing import Optional, Sequence

import scanpy as sc
from anndata import AnnData

import iplot


def rank_marker_genes(adata: AnnData,
                      n_genes: int = 20,
                      groupby: Optional[str] = None,
                      gene_symbols: Optional[str] = None,
                      groups: Optional[Sequence[str]] = None,
                      ret_type: Optional[str] = 'json',
                      save: Optional[str] = None
                      ):
    key = "rank_genes_groups_" + groupby
    data = adata.uns.get(key, None)
    if data is None and groupby is not None:
        sc.tl.rank_genes_groups(adata, groupby=groupby, key_added=key)
        return {'adata': adata,
                "call": f"scanpy.rank_marker_genes(target, groupby=\"{groupby}\", key_added=\"{key}\")",
                "plotly": iplot.scanpy.rank_marker_genes(adata, n_genes=n_genes, gene_symbols=gene_symbols,
                                                         key=key, ret_type=ret_type, save=save, groups=groups)}
    return iplot.scanpy.rank_marker_genes(adata, key=key, n_genes=n_genes, gene_symbols=gene_symbols,
                                          ret_type=ret_type, save=save, groups=groups)
