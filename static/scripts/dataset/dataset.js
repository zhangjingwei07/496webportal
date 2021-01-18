'use strict';

let datasets;
/**
 JS Code for Listing
 **/
const table_option = {
    item: `<tr>
      <td class="name"></td>
      <td class="description"></td>
      <td><span class="n_obs"></span>×<span class="n_vars"></span></td>
      <td class="user"></td>
      <td class="modified"></td>
      <td class="id">
        <button class="btn btn-danger text-white"><i class="fas fa-times"></i></button>
        <button class="btn btn-secondary text-white"><i class="fas fa-cog"></i></button>
      </td>
    </tr>`,
    valueNames: ['name', 'description', 'user', 'n_obs', 'n_vars', 'modified',
        {
            name: 'id',
            attr: 'data-id'
        }
    ],
    page: 25,
    pagination: true
};

const table = new List('dataset-table', table_option, {});

$.get("/dataset/datasets", {
    limit: 10,
    offset: 0
}, data => {
    datasets = data;
    const display = datasets.map(d => {
        d.modified = (new Date(d.modified)).toLocaleString();
        return d;
    });
    table.clear();
    table.add(display);
});

const removeDataset = (id) => {
    if (id === "") {
        return
    }
    $.ajax({
        url: "/dataset/datasets",
        method: "POST",
        data: {
            id: id,
            action: 'DELETE'
        },
    }).done(function (data) {
        table.remove('id', data.id);
    });
};

const updateDataset = () => {
    const id = $(event.target).data("id");
    if (!id) return;
    const data = {id: id, action: "UPDATE"};
    const new_name = $("#new-name").val();
    const new_description = $("#new-description").val();
    if (new_name !== "") {
        data.name = new_name
    }
    if (new_description !== "") {
        data.description = new_description;
    }
    if (!data.name && !data.description) {
        $("#modal-dataset").modal("hide");
        return;
    }
    $.ajax({
        url: "/dataset/datasets",
        method: "POST",
        data: data
    }).done(data => {
        const target = $("td[data-id=" + data.id + "]").parent();
        target.find(".name").text(data.name);
        target.find(".description").text(data.description);
        $("#modal-dataset").modal("hide");
    })
};

const modalDataset = (id) => {
    if (!datasets) {
        return;
    }
    const dataset = datasets.find(dataset => dataset.id === id);
    $("#modal-dataset .name").text(dataset.name);
    $("#modal-dataset .description").text(dataset.description);
    $("#modal-dataset .n_obs").text(dataset.n_obs);
    $("#modal-dataset .n_vars").text(dataset.n_vars);
    $("#modal-dataset .modified").text("last modified: " + dataset.modified.toLocaleString());
    const attrs = JSON.parse(dataset.attrs);
    for (let attr in attrs) {
        const tr_ = $("<tr></tr>");
        const td1_ = $("<td></td>").text(attr);
        const td2_ = $("<td></td>").text(attrs[attr].join("; "));
        tr_.append(td1_, td2_);
        $("#modal-dataset tbody").append(tr_);
    }

    $("#modal-dataset .save").data("id", id);
    $("#modal-dataset").modal('show')
};

$("#dataset-table").click(e => {

    if (!e.target) {
        return;
    }
    let target = $(e.target);
    if (target.hasClass("fas")) {
        target = target.parent();
    }
    if (target.hasClass("btn-danger")) {
        removeDataset(target.parent().data("id"))
    } else if (target.hasClass("btn-secondary")) {
        modalDataset(target.parent().data("id"))
    }
});

$("#modal-dataset").on("hide.bs.modal", () => {
    $("#modal-dataset tr").slice(2).remove();
    $("#new-name").val("");
    $("#new-description").val("")
});


