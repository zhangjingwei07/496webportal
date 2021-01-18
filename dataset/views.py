import importlib
import json
import os
from shutil import rmtree, unpack_archive, get_archive_formats
from time import time

import numpy as np
from django.http import JsonResponse, HttpResponseForbidden, HttpResponseBadRequest
from django.shortcuts import render
from scanpy import read_h5ad

from settings.settings import DATASET_FOLDER, ALLOWED_EXTENSIONS, USER_PROCESS_FOLDER, TEMP_FOLDER
from .models import DataSet
from .utils import get_anndata_attrs


def render_dataset(request):
    return render(request, "dataset/datasets.html")


def render_data_upload(request):
    return render(request, "dataset/dataupload.html",
                  {'allowed_file': ", ".join(ALLOWED_EXTENSIONS)})


def rest_datasets(request):
    # get Datasets
    if request.method == 'GET':
        limit = int(request.GET.get('limit', 0))
        offset = int(request.GET.get('offset', 0))
        result = DataSet.objects.exclude(user__startswith="__")[offset: limit]
        return JsonResponse(list(result.values('id',
                                               'user',
                                               'name',
                                               'description',
                                               'modified',
                                               'n_obs', 'n_vars',
                                               'attrs')), safe=False)
    if request.method == "POST":
        action = request.POST.get('action', "")
        id = request.POST.get("id", None)
        if not id:
            return HttpResponseForbidden()
        # delete datasets
        if action == "DELETE":
            result = DataSet.objects.get(id=id)
            if os.path.isfile(result.path):
                os.remove(result.path)
            if os.path.isdir(result.path):
                rmtree(result.path)
            result.delete()
            return JsonResponse({'id': id})

        # update datasets
        elif action == "UPDATE":
            result = DataSet.objects.get(id=id)
            name = request.POST.get("name")
            if name:
                result.name = name
            description = request.POST.get("description")
            if description:
                result.description = description
            result.save()
            return JsonResponse(result.to_dict(), safe=False)


def data_upload(request):
    """
    Read the uploaded file and then write it as a H5AD file into the file system
    """

    # get file and validate
    file = request.FILES.get('file', None)
    if not file:
        return HttpResponseBadRequest

    file_ori, ext = file.name.rsplit('.', 1)
    ext = ext.lower()
    hash_name = file_ori + "_" + hex(int(time()))[2:]
    path = os.path.join(TEMP_FOLDER, hash_name + "." + ext)

    package = request.POST.get("package")
    method = request.POST.get("method")

    if not method or not package:
        return HttpResponseBadRequest

    with open(path, 'wb+') as f:
        for chunk in file.chunks():
            f.write(chunk)
    # if it is a zipped file, unzip it
    if ext in [x[0] for x in get_archive_formats()]:
        ext_folder = os.path.join(TEMP_FOLDER, file_ori)
        os.mkdir(ext_folder)
        unpack_archive(path, ext_folder)
        os.remove(path)
        path = ext_folder

    try:
        # according to the provided read method, read the data and write as H5AD file
        module = importlib.import_module(package)
        components = method.split(".")
        for attr in components:
            module = getattr(module, attr)
        adata = module(path)
        rmtree(path) if os.path.isdir(path) else os.remove(path)
        data_path = os.path.join(DATASET_FOLDER, hash_name + ".h5ad")
        adata.write(data_path)
    except (FileNotFoundError, ValueError) as e:
        if path:
            rmtree(path) if os.path.isdir(path) else os.remove(path)
        return JsonResponse({'status': False, 'info': 'Invalid data file' + str(e)})
    except Exception as e:
        # if anything wrong, delete anything stored
        if path:
            rmtree(path) if os.path.isdir(path) else os.remove(path)
        return JsonResponse({'status': False, 'info': 'Internal Error: ' + str(e)})

    # update teh file record to database
    saved_file = DataSet(
        user=request.POST.get("owner", "Upload"),
        name=request.POST.get("name", "uploaded-file"),
        path=data_path,
        description=request.POST.get("description", ""),
        n_obs=adata.n_obs,
        n_vars=adata.n_vars,
        attrs=json.dumps(get_anndata_attrs(adata))
    )
    saved_file.save()

    return JsonResponse({'status': True, 'info': "File successfully uploaded as " + saved_file.name + ".h5ad"})


def result_export(request):
    """
    create a new dataset from selected indexes of observations
    """
    pid = request.POST.get("pid", None)
    if not pid:
        return HttpResponseBadRequest
    adata = read_h5ad(os.path.join(USER_PROCESS_FOLDER, str(pid), "results.h5ad"))

    hextime = hex(int(time()))[2:]
    output_path = os.path.join(DATASET_FOLDER, f"exported_{pid}_{hextime}.h5ad")

    indexes = np.fromstring(request.POST.get("index"), dtype=int, sep=",")
    adata = adata[indexes, :]
    adata.write(output_path)

    saved_file = DataSet(
        name=request.POST.get("name", f"export_{pid}"),
        path=output_path,
        description=request.POST.get("description", ""),
        n_obs=adata.n_obs,
        n_vars=adata.n_vars,
        attrs=json.dumps(get_anndata_attrs(adata))
    )
    saved_file.save()
    return JsonResponse(
        {'status': True, 'id': saved_file.id})
