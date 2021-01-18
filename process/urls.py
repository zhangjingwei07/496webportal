from django.urls import path

from . import views

urlpatterns = [
    path('', views.render_new_process),
    path('process.html', views.render_process),
    path('new-process.html', views.render_new_process),
    path('data.html', views.render_data),
    path('process-history.html', views.render_process_history),
    path('process-history', views.get_process_history),
    path('new-process', views.post_new_process)
]
