# Generated by Django 2.2.5 on 2019-12-29 00:52

from django.db import migrations


class Migration(migrations.Migration):

    dependencies = [
        ('process', '0003_auto_20191228_1904'),
    ]

    operations = [
        migrations.RenameField(
            model_name='process',
            old_name='pid',
            new_name='wrid',
        ),
    ]
