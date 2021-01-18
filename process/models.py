from django.db import models


# Create your models here.


class Process(models.Model):
    wrid = models.PositiveSmallIntegerField()
    index = models.PositiveIntegerField()
    call = models.TextField()
    status = models.PositiveSmallIntegerField()
    time = models.DateTimeField(auto_now=True)
    output = models.TextField()
    type = models.CharField(max_length=16)


class WorkerRecord(models.Model):
    """Stores information for the worker"""
    name = models.CharField(max_length=32)
    status = models.PositiveSmallIntegerField()
    curr = models.PositiveSmallIntegerField()
    total = models.PositiveSmallIntegerField()
    time = models.DateTimeField(auto_now=True)
