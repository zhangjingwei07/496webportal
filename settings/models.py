from django.db import models
import json

# Create your models here.
class Methods(models.Model):
    type = models.CharField(max_length=16)
    name = models.CharField(max_length=128)
    package = models.CharField(max_length=64)
    description = models.TextField()
    params = models.TextField()

    def assembly(self):
        return {
            'id': self.id,
            'type': self.type,
            'name': self.name,
            'package': self.package,
            'description': self.description,
            'params': json.loads(self.params)
        }