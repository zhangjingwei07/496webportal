'use strict';

let installedMethods;
const table_option = {
    item: `<tr>
      <td class="type"></td>
      <td class="package"></td>
      <td class="name"></td>
      <td class="description"></td>
      <td>
        <button class="btn btn-danger text-white id delete"><i class="fas fa-times delete"></i></button>
        <button class="btn btn-secondary text-white update"><i class="fas fa-cog update"></i></button>
      </td>
    </tr>`,
    valueNames: ['type', 'package', 'name', 'description', {
        name: 'id',
        attr: 'data-id'
    }
    ],
    page: 10,
    pagination: true
};

const table = new List('methods-table', table_option, {});

$.get("/settings/installed-methods", {
    type: 'processing;reader;plot;iplot',
    name: '_all'
}, data => {
    table.clear();
    installedMethods = data;
    table.add(installedMethods);
});

$("#methods-table tbody").click(() => {

    let target = $(event.target);
    if (target.hasClass("delete")) {
        if (!target.hasClass("id")) {
            target = target.parent();
        }
        $.ajax({
            url: "/settings/installed-methods",
            method: "POST",
            data: {
                id: Number(target.data("id")),
                action: 'DELETE'
            },
        }).done((data) => {
            table.remove('id', data.id);
        });
    } else if (target.hasClass("update")) {
        if (!target.hasClass("btn")) {
            target = target.parent();
        }
        target = target.prev();
        const id = target.data("id");
        const method = installedMethods.find(m => m.id === id);
        $("#type").val(method.type);
        $("#name").val(method.name).removeClass("is-invalid");
        $("#package").val(method.package).removeClass("is-invalid");
        $("#description").val(method.description);
        $("#table-params tbody").empty();
        method.params.forEach(p => insertParamTable(p));
        $("#tab-2").tab('show');
    }
});

$("#div-options").hide();
const parameters = [];
$('#btn-new-param').click(() => {
    $("#modal-new-param").modal("show")
});

$(".btn-checkbox").click(() => {
    $(event.target).toggleClass("btn-success btn-secondary")
});

$("#param-type").on('change', () => {
    if ($(event.target).val() === "text") {
        $("#param-is-list").show();
    } else {
        $("#param-is-list").hide();
    }
    if ($(event.target).val() === "option") {
        $("#div-options").show();
    } else {
        $("#div-options").hide();
    }
});

$(".required").on("change", () => {
    if ($(event.target).val() !== "") {
        $(event.target).removeClass("is-invalid")
    } else {
        $(event.target).addClass("is-invalid")
    }
});

const addParam = () => {
    const param = {
        name: $("#param-name").val(),
        type: $('#param-type :selected').text()
    };
    if (param.name === "") {
        return
    }
    const annotation = $("#param-annotation").val().trim();
    if (annotation !== "") {
        param.annotation = annotation;
    }
    if (param.type === "option") {
        const options = $('#param-options').val().trim();
        if (options === "") {
            return
        }
        param.options = options.split(",")
    } else if (param.type === "text" && $("#param-is-list").hasClass("btn-success")) {
        param.isList = true;
    }
    if ($("#param-required").hasClass("btn-success")) {
        param.required = true;
    }
    const default_div = $("#param-default");
    if (default_div.hasClass("is-invalid")) {
        return;
    }
    if (default_div.val() !== "") {
        param.default = getDefaultValue(default_div);
    }
    insertParamTable(param);
    parameters.push(param);
    $("#form-param").trigger("reset");
    $("#param-name").addClass("is-invalid");
    $("#modal-new-param").modal("hide");
};

const getDefaultValue = (default_div) => {
    $("#param-default").removeClass("is-invalid");
    let v = default_div.val().trim();
    const type = $('#param-type :selected').text();
    if (type === "number") {
        v = Number(v);
    } else if (type === "option" && ($("#param-options").val().split(",").includes(v) === false)) {
        console.log("www");
        console.log(v);
        v = NaN;
    } else if (type === "bool") {
        if (v.toLowerCase() === "true") {
            v = true;
        } else if (v.toLowerCase() === "false") {
            v = false;
        } else {
            v = NaN;
        }
    }
    if (isNaN(v) && typeof v !== "string") {
        $("#param-default").addClass("is-invalid");
        return;
    }
    if (v === "EMPTY_STRING") {
        v = ""
    }
    return v;
};
$("#param-default").on("change", () => {
    getDefaultValue($("#param-default"))
});

const insertParamTable = (param) => {
    const _tr = $("<tr></tr>");

    _tr.append($("<td class='name'></td>").text(param.name));
    _tr.append($("<td class='type'></td>").text(param.type));
    if (!param.annotation) {
        _tr.append($("<td class='annotation'></td>"))
    } else {
        _tr.append($("<td class='annotation'></td>").text(param.annotation))
    }
    _tr.append($('' +
        '<td>' +
        '<button type="button" class="btn btn-danger my-1 text-white delete">' +
        '<i class="fas fa-times delete"></i>' +
        '</button>' +
        '</td>'));
    $("#table-params tbody").append(_tr)
};

$("#table-params tbody").click(() => {
    if (!$(event.target).hasClass("delete")) {
        return;
    }
    let target = $(event.target).parent().parent();
    if (target.prop("tagName") !== "TR") {
        target = target.parent();
    }
    parameters.splice(parameters.findIndex((e) => {
        return e.name === target.find(".name").text()
    }), 1);
    target.remove();
});

const postMethod = () => {
    const method_data = {
        type: $("#type").val(),
        package: $("#package").val(),
        name: $("#name").val(),
        description: $("#description").val(),
        params: JSON.stringify(parameters)
    };
    if (method_data.package !== "" && method_data.name !== "") {
        $.ajax({
            url: '/settings/update-method',
            data: method_data,
            type: 'POST',
            success: data => {
                window.location.reload();
            }
        });
    }
};

const resetForm = () => {
    $("#form-method").trigger("reset");
    parameters.splice(0, parameters.length);
    $("#table-params tbody").empty();
};

const exportMethods = () => {
    const data = {
        name: "SCAWPMETHODS"
    };
    data.data = installedMethods.map(m => {
        const new_m = m;
        delete new_m.id;
        return new_m;
    });
    download(JSON.stringify(data), "installed_methods.json");
};

const importMethods = () => {
    const formData = new FormData($("#form-upload")[0]);
    $.ajax({
        url: '/settings/reset-methods',
        data: formData,
        processData: false,
        contentType: false,
        type: 'POST',
        success: data => {
            window.location.reload();
        }
    })
};

