from anndata import AnnData


def get_anndata_attrs(adata: AnnData):
    attrs = {}
    for attr in ["obs", "var", "uns", "obsm", "varm", "layers"]:
        keys = getattr(adata, attr).keys()
        if len(keys) > 0:
            attrs[attr] = list(keys)
    return attrs
