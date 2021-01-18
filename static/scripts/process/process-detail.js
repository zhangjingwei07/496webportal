'use strict';

const report_div = $("#report-div");
const url = {
    name: new URLSearchParams(window.location.search).get("id")
};
$.get("/process-history", url, (data) => {
    data.sort((a, b) => {
        return -1 * (b.index - a.index);
    });
    data.forEach((p) => {
        const col_ = $('<div class="col-md-4 col-lg-3"></div>');
        const card_ = $('<div class="card-progress-report card-border mb-3 card card-body"></div>');
        if (p.status === 1) {
            card_.addClass("border-success");
        } else if (p.status === 0) {
            card_.addClass("border-primary");
        } else {
            card_.addClass("border-danger");
        }
        const head_ = $('<h5 class="text-lowercase card-title"></h5>');
        if (p.index === 0) {
            head_.text("target = " + p.call.replace("target", "filename"))
        } else {
            head_.text(p.call)
        }
        const time_ = $('<h6 class="card-subtitle"></h6>');
        const output_ = $("<p></p>");
        output_.append(p.output.replace(/\n/g, "<br>"));
        const time = new Date(p.time);
        time_.text("Updated on: " + time.toLocaleString());
        card_.append(head_, time_, output_);
        col_.append(card_);
        report_div.append(col_)
    })
});

$("#script-button").click(() => {
    let script_ = "";
    $("#report-div .card-title").each((index, obj) => {
        script_ += ($(obj).text() + "\n");
    });
    download(script_, $("#worker-name").text() + "_script.py", "text/plain");
});

$("#report-button").click(() => {
    let report_ = "";
    const cards = $("#report-div .card");
    cards.each((index, obj) => {
        report_ += $(obj).find(".card-title").text() + "\n";
        report_ += $(obj).find(".card-subtitle").text() + "\n";
        report_ += $(obj).find("p").text() + "\n\n"
    });
    download(report_, $("#worker-name").text() + ".txt", "text/plain");
});


const remove_history = (id) => {
    $.ajax({
        url: '/process/process-history',
        data: {id: id, action: 'DELETE'},
        type: 'POST',
        success: (data) => {
            window.location.href = "/process/process-history.html"
        }
    });
};
