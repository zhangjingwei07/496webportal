from django.urls import path

from . import views

urlpatterns = [
    path('', views.render_installed_methods),
    path('installed-methods.html', views.render_installed_methods),
    path('installed-methods', views.get_installed_methods),
    path('update-method', views.update_installed_methods),
    path('reset-methods', views.reset_methods)
]
