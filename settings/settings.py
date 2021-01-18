import os
from shutil import get_archive_formats

from SIMLR.settings import BASE_DIR

USER_DATA_DIR = os.path.join(BASE_DIR, "userData")
if not os.path.exists(USER_DATA_DIR):
    os.mkdir(USER_DATA_DIR)

ALLOWED_EXTENSIONS = set(['h5ad', 'csv', 'h5', 'loom', 'mtx', 'txt'] + [x[0] for x in get_archive_formats()])
DATASET_FOLDER = os.path.join(USER_DATA_DIR, "datafile")
if not os.path.exists(DATASET_FOLDER):
    os.mkdir(DATASET_FOLDER)

TEMP_FOLDER = os.path.join(USER_DATA_DIR, "temp")
if not os.path.exists(TEMP_FOLDER):
    os.mkdir(TEMP_FOLDER)

USER_PROCESS_FOLDER = os.path.join(USER_DATA_DIR, "processes")
if not os.path.exists(USER_PROCESS_FOLDER):
    os.mkdir(USER_PROCESS_FOLDER)
