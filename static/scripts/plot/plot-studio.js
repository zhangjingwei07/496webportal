'use strict';
const shape = ($("#dataset-shape").text()).split(", ").map(e => Number(e));
const n_vars = shape[0];
const n_obs = shape[1];
const attrs = JSON.parse($("#dataset-attrs").text());
const id = Number($(".id").text());

const cluster_methods = ['louvain', 'leiden'];


const resetStudio = () => {
    $("#plotly").prop("hidden", true);
    $("#div-select").prop("hidden", false);
    deactivateOptions();
};

if (attrs.obsm.includes("X_pca")) {
    $(".basis").append($("<option value='iplot.scanpy.pca'>PCA</option>"))
}
if (attrs.obsm.includes("X_tsne")) {
    $(".basis").append($("<option value='iplot.scanpy.tsne'>tSNE</option>"))
}
if (attrs.obsm.includes("X_umap")) {
    $(".basis").append($("<option value='iplot.scanpy.umap'>UMAP</option>"))
}

attrs.obs.filter(col => cluster_methods.includes(col)).forEach(col => {
    $(".groups").append($("<option>" + col + "</option>"))
});


$("#scatter-plot").click(() => {
    const call = $("#scatter-basis").val();
    const names = $("#scatter-names").selectpicker("val");
    if (!call || !names || names.length === 0) {
        return;
    }
    $.post("/plot/plot-sync", {
        id: id,
        call: call,
        params: JSON.stringify({names: names})
    }, data => {
        $("#plotly").prop("hidden", false);
        Plotly.newPlot('plotly', JSON.parse(data), {});
        activateOptions("scatter");
        $("#div-select").prop("hidden", true);
        window.scrollTo(0, document.body.scrollHeight);
    });
});

$("#rank-marker-genes-plot").click(() =>{
    const call = $("#rank-marker-genes-method").val();
    const groupby = $("#rank-marker-genes-names").val();
    $.post("/plot/plot-sync", {
        id: id,
        call: call,
        params: JSON.stringify({groupby: groupby})
    }, data => {
        $("#plotly").prop("hidden", false);
        Plotly.newPlot('plotly', JSON.parse(data), {});
        activateOptions("rank-marker-genes");
        $("#div-select").prop("hidden", true);
        window.scrollTo(0, document.body.scrollHeight);
    });
});
