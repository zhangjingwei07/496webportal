import json
import os
from shutil import rmtree

from django.http import JsonResponse, HttpResponseBadRequest
from django.shortcuts import render, get_object_or_404
from scanpy import read_h5ad

from dataset.models import DataSet
from settings.settings import USER_PROCESS_FOLDER
from .models import *
from .worker import Worker


def render_new_process(request):
    return render(request, "process/new-process.html")


def render_process(request):
    id_ = int(request.GET.get('id', None))
    if id_ is None:
        return HttpResponseBadRequest
    worker = get_object_or_404(WorkerRecord, id=int(id_))
    return render(request, "process/process.html", {'worker': worker})


def render_process_history(request):
    return render(request, "process/process-history.html")


def render_data(request):
    id_ = int(request.GET.get('id', None))
    if id_ is None:
        return HttpResponseBadRequest
    worker = get_object_or_404(WorkerRecord, id=int(id_))
    dataset = get_object_or_404(DataSet, name=f"Worker_{id_}")
    dataset.attrs = json.loads(dataset.attrs)
    path = os.path.join(USER_PROCESS_FOLDER, str(worker.id), 'results.h5ad')
    annData = read_h5ad(path)
    return render(request, "process/data.html",
                  {'worker': worker,
                   'vars': annData.var.reset_index().to_html(index=False,
                                                             classes='mb-0 table table-bordered',
                                                             max_rows=10000),
                   'obs': annData.obs.reset_index().to_html(index=False,
                                                            classes='mb-0 table table-bordered',
                                                            max_rows=10000),
                   'dataset': dataset
                   }
                  )


def get_process_history(request):
    if request.method == 'GET':
        name = request.GET.get("name", "")
        if name == "_all":
            result = WorkerRecord.objects.all()
        else:
            result = Process.objects.filter(wrid=name)
        return JsonResponse(list(result.values()), safe=False)

    if request.method == 'POST' and request.POST.get('action') == 'DELETE':
        id = request.POST.get("id", "")
        WorkerRecord.objects.filter(id=id).delete()
        Process.objects.filter(wrid=id).delete()
        DataSet.objects.filter(name=f'Worker_{id}').delete()
        try:
            rmtree(os.path.join(USER_PROCESS_FOLDER, id))
        except FileNotFoundError:
            pass
        return JsonResponse({'id': id, 'status': 'success'})


def post_new_process(request):
    process = json.loads(request.POST.get('process'))
    worker = Worker(process, request.POST.get('name'))
    integrity = worker.check_integrity()
    if integrity['status']:
        worker.start()
    return JsonResponse(integrity)
