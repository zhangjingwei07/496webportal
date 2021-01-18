const allowed_file = $("#allowed_file").text().split(", ");
let installed_readers = [];
$.get("/settings/installed-methods", {
    type: 'reader',
    name: '_all'
}, data => {
    const options = $('#reader');
    data.forEach(reader => {
        const option = $('<option>' + reader.package + '.' + reader.name + '</option>');
        option.data('id', reader.id);
        options.append(option);
        installed_readers.push(reader);
    });

});

/***
 * Update the file information when user choose data file
 */
$("#file").change(() => {
    if ($("#file")[0].files.length === 0) {
        return;
    }
    const file = $("#file")[0].files[0];
    if (file.name.lastIndexOf(".") === -1) {
        return;
    }
    const file_name = file.name.substring(0, file.name.lastIndexOf("."));
    const file_ext = file.name.substring(file.name.lastIndexOf(".") + 1, file.name.length);

    $("#file-text").text("File size: " + String(Math.round10(file.size / 1024 / 1024, -2)) + "MB");
    $("#name").val(file_name);
});


/***
 * Reset the Whole Form
 */
const resetForm = () => {
    $('#form-upload').trigger('reset');
    $("#file-text").text("For .mtx, please upload the directory as a zip file")
};

/***
 * Make the Ajax call that upload the file to the server
 */
const uploadFile = () => {
    const file = $("#file")[0].files[0];
    if (!file) {
        return;
    }
    const file_ext = file.name.substring(file.name.lastIndexOf(".") + 1, file.name.length);
    if (!allowed_file.includes(file_ext)) {
        $("#modal-warning .modal-title").text("Incompatible file Extension");
        $("#modal-warning .modal-body p").text("The uploading file does not have a recognizable extension, " +
            "please upload a valid dataset file");
        $("#modal-warning").modal();
        return
    }
    const reader = installed_readers.find(e => {
        return e.id === $('#reader option:selected').data('id')
    });

    const formData = new FormData($("#form-upload")[0]);
    formData.append('name', $('#name').val());
    formData.append('owner', $('#owner').val());
    formData.append('description', $('#description').val());
    formData.append('package', reader.package);
    formData.append('method', reader.name);

    $("#modal-warning .modal-title").text("Uploading");
    $("#modal-warning .modal-body p").text("The file is uploading, please keep this page open");
    $("#modal-warning").modal();

    $.ajax({
        url: '/dataset/data-upload',
        data: formData,
        processData: false,
        contentType: false,
        type: 'POST',
        success: data => {
            $("#modal-warning .modal-title").text("Uploaded");
            $("#modal-warning .modal-body p").text(data.info);
        }
    });
};