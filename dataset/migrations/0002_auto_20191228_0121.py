# Generated by Django 2.2.5 on 2019-12-28 06:21

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('dataset', '0001_initial'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='datafile',
            name='path',
        ),
        migrations.AddField(
            model_name='datafile',
            name='file',
            field=models.FileField(null=True, upload_to='D:\\Repo\\SIMLR\\datafile'),
        ),
    ]
