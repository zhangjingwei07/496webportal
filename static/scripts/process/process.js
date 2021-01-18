'use strict';
/**
 For Process Section
 **/
$("#submit-footer, #active-process-table, #new-steps-card").hide();

/** ===== Table Options ===== **/
const new_steps_option = {
    item: `<div class="col-lg-2 col-md-3">
            <div class="card-shadow-secondary border mb-3 card card-body border-secondary card-new-steps">
              <h5 class="card-title" style="text-transform: lowercase;">
                <span class="package"></span>.<span class="name">
              </h5>
              <div class="scroll-area-sm">
                <div class="scrollbar-container ps--active-y">
                  <span class="description"></span>
                </div>
              </div>
            </div>
          </div>`,
    valueNames: ['name', 'package', 'description'],
    page: 12,
    pagination: true
};


const dataset_table_option = {
    item: `<tr>
    <td class="name id"></td>
    <td><span class="n_obs"></span>×<span class="n_vars"></span></td>
    <td class="description"></td>
    <td class="modified">
    </td>
  </tr>`,
    valueNames: ['name', 'description', 'modified', "n_obs", "n_vars",
        {
            name: 'id',
            attr: 'data-id'
        }
    ],
    page: 10,
    pagination: true
};

/** ===== Variables ===== **/
const new_steps_table = new List('new-steps-table', new_steps_option, {});

const active_processing = [];

let installedMethods;

$.get("/settings/installed-methods", {
    type: 'processing;plot;iplot',
    name: '_all',
}, data => {
    installedMethods = data;
    new_steps_table.clear();
    new_steps_table.add(installedMethods);
});

$("#new-steps-table .list").click(() => {
    let target = $(event.target);
    for (let i = 0; i < 5; i++) {
        if (!target.hasClass("card-new-steps")) {
            target = target.parent();
        } else {
            break;
        }
    }
    if (!target.hasClass("card-new-steps")) {
        return;
    }
    const process_info = constructProcess(target.find(".name").text(),
        target.find(".package").text());
    active_processing.push(process_info);
});


const dataset_table = new List('dataset-table', dataset_table_option, {});

$.get("/dataset/datasets", {
    limit: 10,
    offset: 0
}, data => {
    dataset_table.clear();
    dataset_table.add(data);

    const url = new URL(window.location.href);
    const dataset = Number(url.searchParams.get("dataset"));
    if (dataset) {
        const found_dataset = dataset_table.items.find(data => data._values.id === dataset);
        if (found_dataset) {
            const chosen_dataset = $("#chosen-dataset");
            chosen_dataset.val(found_dataset._values.name);
            chosen_dataset.data("id", found_dataset._values.id);
            $("#modal-dataset").modal("hide");

            const input_ = $('#chosen-dataset');
            input_.removeClass("is-invalid");
            $("#submit-footer, #active-process-table, #new-steps-card").show();
        }
    }
});


$("#dataset-table .list").click(() => {
    let target = $(event.target);
    if (target.prop("tagName") === "TD") {
        target = target.parent()
    } else {
        return
    }
    const chosen_dataset = $("#chosen-dataset");
    chosen_dataset.val($(target.children()[0]).text());
    chosen_dataset.data("id", $(target.children()[0]).data("id"));
    $("#modal-dataset").modal("hide");

    const input_ = $('#chosen-dataset');
    input_.removeClass("is-invalid");
    $("#submit-footer, #active-process-table, #new-steps-card").show();
});


$("#active-process-table").click(e => {
    let target = $(event.target);
    if (display_mode === "s"
        && !target.hasClass("fas")
        && !target.hasClass("option-btn")) {
        if (target.parent().hasClass("card-pp")) {
            target = target.parent();
        }
        if (!target.hasClass("card-pp")) {
            return;
        }
        target.children("span").toggle();
        target.parent().toggleClass("col-lg-2 col-md-3 col-sm-6 col-lg-1 col-md-2 col-sm-4");
        return;
    }
    if (target.hasClass("fas")) {
        target = target.parent()
    }
    const pid = target.parent().data("pid");
    if (pid !== 0 && !pid) {
        return
    }
    if (target.hasClass("btn-alternate")) {
        const params = active_processing.find(el => el.pid === pid).params;
        const modal_option = $("#modal-option");
        const option_content = $("#option-content");
        const option_bool = $("#option-bool");
        modal_option.data("type", "processing");
        option_content.empty();
        option_bool.empty();
        modal_option.data("pid", target.parent().data("pid"));
        constructOptions(option_content, option_bool, params);
        modal_option.modal("show");
    } else if (target.hasClass("btn-danger")) {
        const index_ = active_processing.findIndex(el => el.pid === pid);
        if (index_ !== -1) {
            active_processing.splice(index_, 1)
        }
        for (let i = 0; i < 5; i++) {
            if (target.hasClass("sort")) {
                target.remove();
                return;
            } else {
                target = target.parent()
            }
        }
    } else if (target.hasClass("btn-secondary")) {
        active_processing.find(el => el.pid === pid).view = true;
        target.removeClass("btn-secondary");
        target.addClass("btn-success")
    } else if (target.hasClass("btn-success")) {
        active_processing.find(el => el.pid === pid).view = false;
        target.removeClass("btn-success");
        target.addClass("btn-secondary")
    }
});


$("#modal-option").on("hide.bs.modal", () => {
    const target = $("#modal-option");
    let obj;
    if (target.data("type") === "processing") {
        obj = active_processing.find(el => el.pid === target.data("pid"))
    }
    if (!obj) {
        return
    }
    let required_params = false;
    for (let i = 0; i < obj.params.length; i++) {
        obj.params[i].default = retriveValue(obj.params[i]);
        if (obj.params[i].required && obj.params[i].default === "") {
            required_params = true;
        }
    }
    const card_ = $("#active-process-table .card-pp span")
        .filter((index, e) => $(e).data("pid") === obj.pid).parent();

    if (required_params) {
        card_.addClass("card-shadow-danger border-danger");
        card_.removeClass("card-shadow-secondary border-secondary");
    } else {
        card_.removeClass("card-shadow-danger border-danger");
        card_.addClass("card-shadow-secondary border-secondary");
    }

});

/** hide the processes details in the queue
 **/
let display_mode = "l";

$("#change-queue-display").click(() => {
    $("#change-queue-display").toggleClass("btn-primary btn-secondary");
    const targeted_cards = $("#active-process-table .card-pp");
    if (display_mode === "l") {
        display_mode = "s";
        targeted_cards.parent()
            .removeClass("col-lg-2 col-md-3 col-sm-6");
        targeted_cards.parent()
            .addClass("col-lg-1 col-md-2 col-sm-4");
        $("#active-process-table .scroll-area-sm.mt-1").hide();
        $("#active-process-table .card-pp span").hide();
        $("#active-process-table .card-process").css("height", "128px");
    } else {
        display_mode = "l";
        targeted_cards.parent()
            .addClass("col-lg-2 col-md-3 col-sm-6");
        targeted_cards.parent()
            .removeClass("col-lg-1 col-md-2 col-sm-4");
        $("#active-process-table .scroll-area-sm.mt-1").show();
        $("#active-process-table .card-pp span").show();
        $("#active-process-table .card-process").css("height", "240px");
    }
});
