from django.urls import path

from . import views

urlpatterns = [
    path('', views.render_plots),
    path('plots.html', views.render_plots),
    path('plot-detail.html', views.render_plot_detail),
    path('plot-studio.html', views.render_plot_studio),
    path('plot-sync', views.render_plot_studio_detail)
]
