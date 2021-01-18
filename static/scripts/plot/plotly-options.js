let selected_points = [];
$(".scatter").hide();
const divPlotly = $("#plotly")[0];

const activateOptions = (type) => {
    $("#orientation").click(() => {
        const target = $(event.target);
        if (target.data('orientation') === 'v') {
            Plotly.restyle(divPlotly, {orientation: 'h'});
            target.data('orientation', 'h').text("Vertical Layout")
        } else {
            Plotly.restyle(divPlotly, {orientation: 'v'});
            target.data('orientation', 'v').text("Horizontal Layout")
        }
    });

    $("#legend").click(() => {
        const legend = $("#legend");
        if (legend.text() === "Hide Legend") {
            legend.html('Show Legend');
            Plotly.relayout(divPlotly, {showlegend: false});
        } else {
            legend.html('Hide Legend');
            Plotly.relayout(divPlotly, {showlegend: true});
        }
    });

    const templateOption = $("#template");
    templateOption.click(() => {
        let template = templateOption.data('template');
        if (!template) {
            templateOption.data("template", "plotly_white");
            template = "plotly_white";
        }
        updateTemplate(template);
        if (template === "plotly_white") {
            templateOption.data('template', 'plotly_dark').text("Dark Background");
        } else {
            templateOption.data('template', 'plotly_white').text("Bright Background");
        }
    });

    const updateTemplate = (template) => {
        $.ajax({
            url: '/static/plot-template/' + template + '.json',
            dataType: 'json',
            success: (data) => {
                Plotly.relayout(divPlotly, {template: data});
            },
            statusCode: {
                404: () => alert('Plot not found')
            }
        });
    };

    if (type === "scatter") {
        $(".scatter").show();
        divPlotly.on('plotly_selected', (data) => {
            selected_points = [];
            if (data) {
                selected_points.push(...(data.points.map(p => p.text)));
            }
        });
        $("#select-points").click(() => {
            $("#modal-data-wizard").data("type", "points").modal("show");
            $("#selected-length").text(selected_points.length);
        });
        $("#select-traces").click(() => {
            $("#modal-data-wizard").data("type", "traces").modal("show");
            const vis = divPlotly.data.filter(trace => trace.visible === true || trace.visible === undefined);
            let sum = 0;
            vis.forEach(trace => sum = sum + trace.text.length);
            $("#selected-length").text(sum);
        });
        $("#select-confirm").prop("disabled", false).click(() => {
            const modal_data_wizard = $("#modal-data-wizard");
            const type = modal_data_wizard.data("type");
            let indexes;
            if (type === "points") {
                indexes = selected_points;
            } else if (type === "traces") {
                indexes = [];
                const vis = divPlotly.data.filter(trace => trace.visible === true || trace.visible === undefined);
                vis.forEach(trace => indexes.push(...trace.text));
            } else {
                return;
            }
            const data = {
                pid: $(".id").text(),
                index: indexes.join(",")
            };

            const name = $("#name-input").val();
            if (name !== "") data.name = name;
            const description = $("#description-input").val();
            if (description !== "") data.description = description;

            $("#modal-warning .modal-title").text("Exporting");
            $("#modal-warning .modal-body p").text("The data is exporting, please keep this page open");
            $("#modal-warning").modal();
            modal_data_wizard.modal("hide");
            $.post('/dataset/result-export', data, (data) => {
                    $("#modal-warning .modal-title").text("Export");
                    $("#modal-warning .modal-body p").text("Success");
                    const href = "/process/new-process.html?dataset=" + String(data.id);
                    const further_button = $("<a class='btn btn-primary'>Start Further Analysis</a>");
                    further_button.attr("href", href);
                    $("#modal-warning .modal-footer").prepend(further_button);
                }
            );
        });
        $("#rank_marker_genes").click(() => {
            const data = {
                groupby: divPlotly.layout.title.text.split(" - "),
                groups: divPlotly.data.filter(trace => trace.visible === true).map(trace => trace.name)
            }
            if (data.groupby.length === 2) {
                data.groupby = data.groupby[1]
            } else {
                return;
            }
            $.post("/plot/plot-sync", {
                id: id,
                call: "studio.rank_marker_genes",
                params: JSON.stringify(data)
            }, data => {
                Plotly.newPlot('plotly-sub', JSON.parse(data), {});
                $("#modal-plotly").modal("show");
            });
        })
    }
};


const deactivateOptions = (type) => {
    $(".plotly-options").off();
    $("#select-confirm").off().prop("disabled", true);
    selected_points = [];
    $(".scatter").hide();
};