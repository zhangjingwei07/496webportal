const divPlotly = document.getElementById('plotly');
let plot;
let selected;
let selected_points = [];
const templates = {
    none: {}
};
$(document).ready(() => {
    $.ajax({
        url: $('#output').data('link'),
        dataType: 'json',
        success: function (data) {
            Plotly.newPlot('plotly', JSON.parse(data), {});
            activateOptions(divPlotly);
            divPlotly.on('plotly_selected', (data) => {
                selected = data;
            });
        },
        statusCode: {
            404: () => alert('Plot not found')
        }
    });


});