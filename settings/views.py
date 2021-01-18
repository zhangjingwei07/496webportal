import json

from django.http import JsonResponse, Http404, HttpResponseForbidden, HttpResponseBadRequest
from django.shortcuts import render

from .models import Methods


def render_installed_methods(request):
    return render(request, "settings/installed-methods.html")


def get_installed_methods(request):
    if request.method == "GET":
        type = request.GET.get("type", "").split(";")
        name = request.GET.get("name", "")
        if name == "_all":
            read_json = Methods.objects.filter(type__in=type)
        else:
            read_json = Methods.objects.filter(name=name, type__in=type)
        read_json = [m.assembly() for m in read_json]
        return JsonResponse(read_json, safe=False)
    if request.method == "POST":
        action = request.POST.get('action', "")
        id = request.POST.get("id", None)
        if action != "DELETE" or not id:
            return HttpResponseForbidden()
        Methods.objects.get(id=id).delete()
        return JsonResponse({'id': id})
    return Http404


def update_installed_methods(request):
    saved_method, _ = Methods.objects.update_or_create(
        type=request.POST.get('type'),
        package=request.POST.get('package'),
        name=request.POST.get('name'),
        defaults={
            'description': request.POST.get('description'),
            'params': str(request.POST.get('params', ""))
        }
    )
    return JsonResponse(saved_method.assembly(), safe=False)


def reset_methods(request):
    file = request.FILES.get('file', None)
    if not file:
        return HttpResponseBadRequest
    data = json.load(file)
    if data.get("name", "") != "SCAWPMETHODS":
        return HttpResponseBadRequest
    methods = data.get("data", "")
    Methods.objects.all().delete()
    bulk = []
    for method in methods:
        bulk.append(Methods(
            type=method['type'],
            name=method['name'],
            package=method['package'],
            description=method.get("description", ""),
            params=json.dumps(method.get('params'))
        ))
    Methods.objects.bulk_create(bulk)
    return JsonResponse({'info': 'imported'})
