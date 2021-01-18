from django.db import models


class DataSet(models.Model):
    user = models.CharField(max_length=32)
    name = models.CharField(max_length=256)
    description = models.TextField()
    path = models.FilePathField(null=True)
    modified = models.DateTimeField(auto_now=True)
    n_obs = models.IntegerField()
    n_vars = models.IntegerField()
    attrs = models.TextField()

    def to_dict(self):
        return {
            'id': self.id,
            'user': self.user,
            'name': self.name,
            'description': self.description,
            'modified': self.modified,
            'n_obs': self.n_obs,
            'n_vars': self.n_vars,
            'attrs': self.attrs
        }
