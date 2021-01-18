# Single Cell Analysis Web Portal

CSC494 Independent Research Project.  
Under the supervision of Prof. Bo Wang

## Demo
This video demo shows the basic workflow 

https://www.youtube.com/watch?v=bbk_lKrU3No&t=12s

## Install

```
conda env create -f environment.yml
python manage.py runserver
```

IMPORTANT: There is a dependency issue among h5py, annData packages. DON'T upgrade any packages yet. 

Major dependencies: django, scanpy, plotly, plotly-orca, requests

## Contribute
See the development decumentation on https://github.com/lihd1003/Single-Cell-Analysis-Web-Portal/blob/master/docs/README.md
