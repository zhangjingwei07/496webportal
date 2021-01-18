$("#submit-process").click(e => {
    const worker_name = $("#worker-name");
    if (worker_name.val() === "") {
        return;
    }
    const full_data = {
        name: worker_name.val()
    };
    const process_order = $("#active-process-table").find(".btn-alternate.option-btn")
        .parent()
        .map((index, e) => {
            return $(e).data("pid")
        })
        .get();

    const reader_data = {
        name: "read_h5ad",
        package: "scanpy",
        params: {filename: $("#chosen-dataset").data("id")},
        type: 'reader'
    };

    if (!reader_data.params.filename) {
        $("#chosen-dataset").addClass('is-invalid');
        return;
    }

    const data = [reader_data];

    let integrity = true;
    process_order.forEach((order) => {
        const process = active_processing.find(e => e.pid === order);
        const target = {
            name: process.name,
            type: process.type,
            package: process.package,
            params: {}
        };
        if (process.view) {
            target['view'] = true
        }

        process.params.forEach(p => {
            if (p.required && p.default === "") {
                integrity = false;
                const card_ = $("#active-process-table .card-pp span")
                    .filter((index, e) => $(e).data("pid") === order).parent();
                card_.addClass("card-shadow-danger border-danger");
                card_.removeClass("card-shadow-secondary border-secondary");
            }
            target.params[p.name] = p.default;
            if (p.type === "number") {
                target.params[p.name] = Number(target.params[p.name])
            }
            if (p.isList && p.default) {
                target.params[p.name] = p.default.split(", ");
            }
        });
        data.push(target);
    });

    if (!integrity) {
        $("#modal-warning .modal-title").text("Integrity Check");
        $("#modal-warning .modal-body p").text("Exist required parameters that is not filled");
        $("#modal-warning").modal();
        return;
    }

    full_data.process = JSON.stringify(data);

    console.log(full_data);

    $.ajax({
        url: '/process/new-process',
        data: full_data,
        type: 'POST',
        success: data => {
            $("#modal-warning .modal-title").text("Work Deplyed");
            $("#modal-warning .modal-body p").text("Work has been successfully deployed, " +
                "the deployment ID is " + data.info);
            $("#modal-warning").modal();
        }
    });
});
