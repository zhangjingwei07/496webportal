from django.urls import path

from . import views

urlpatterns = [
    path('dataupload.html', views.render_data_upload),
    path('datasets.html', views.render_dataset),
    path('datasets', views.rest_datasets),
    path('data-upload', views.data_upload),
    path('result-export', views.result_export)
]
