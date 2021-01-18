import importlib
import json
import os
import traceback
import torch
from torchvision import models
from torch import nn
import pytorch_lightning as pl
from sklearn.model_selection import train_test_split
from torch.utils.data import DataLoader, random_split, Subset, Dataset
from torch.nn import functional as F
from torchvision import datasets, transforms
import torch.optim as optim
from tqdm import tqdm
from pytorch_lightning.callbacks import ModelCheckpoint
from collections import Counter
import statistics

from torchvision.datasets.utils import download_and_extract_archive
import pandas as pd
import numpy as np

from settings.settings import USER_PROCESS_FOLDER
from .models import Process

###------------------------------Model---------------------------------------###
TOP_K = 5

class CSQLightening(pl.LightningModule):
  def __init__(self,n_class,n_features,batch_size=64,l_r=1e-5,lamb_da=0.0001,beta=0.9999,bit=64):
    super(CSQLightening, self).__init__()
    print("hparam: l_r = {}, lambda = {}, beta = {}", l_r, lamb_da, beta)
    self.batch_size = batch_size
    self.l_r = l_r
    self.bit = bit
    self.n_class = n_class
    self.lamb_da = lamb_da
    self.beta = beta
    self.samples_in_each_class = None # Later initialized in training step
    self.hash_centers = get_hash_centers(self.n_class, self.bit)
    ##### model structure ####
    # input size = batch size * 19791
    self.hash_layer = nn.Sequential(
        nn.Linear(n_features, 6300),
        nn.ReLU(inplace=True),
        nn.Dropout(0.2),
        nn.Linear(6300, 2100),
        nn.ReLU(inplace=True),
        nn.Dropout(0.2),
        nn.Linear(2100, 710),
        nn.ReLU(inplace=True),
        nn.Dropout(0.2),
        nn.Linear(710, 200),
        nn.ReLU(inplace=True),
        nn.Linear(200, self.bit),
    )
    

  def forward(self, x):
    # forward pass returns prediction
      x = self.hash_layer(x)
      return x

  def CSQ_loss_function(self, hash_codes, labels):
    hash_codes = hash_codes.tanh()
    hash_centers = self.hash_centers[labels]
    hash_centers = hash_centers.type_as(hash_codes)

    # Class-Balanced Loss on Effective Number of Samples
    # Reference Paper https://arxiv.org/abs/1901.05555
    if self.samples_in_each_class == None:
      self.samples_in_each_class = self.trainer.datamodule.samples_in_each_class
      self.n_class = self.trainer.datamodule.N_CLASS
    class_sample_count = self.samples_in_each_class[labels]
    weight = (1 - self.beta)/(1 - torch.pow(self.beta, class_sample_count))
    weight = weight / weight.sum() * self.n_class
    weight = weight.type_as(hash_codes)

    # Center Similarity Loss
    BCELoss = nn.BCELoss(weight=weight.unsqueeze(1).repeat(1,self.bit))
    C_loss = BCELoss(0.5 * (hash_codes + 1),
                        0.5 * (hash_centers + 1))
    # Quantization Loss
    Q_loss = (hash_codes.abs() - 1).pow(2).mean()

    loss = C_loss + self.lamb_da * Q_loss
    return loss



class WorkerStep:
    """
    One step in a worker, manage the logging and running of one step
    @:param context: an dictionary object contains the function call
    @:param call: a string of the function call
    @:param output: a string of the function call's output
    @:param status: the status of this step
    @:param file: the associated file of this step, if needed
    """

    def __init__(self, context, wrid, index, file, annData):
        self.context = context
        self.output = ""
        self.status = 0
        self.wrID = wrid
        self.index = index
        self.file = file
        self.annData = annData
        self.folder = os.path.join(USER_PROCESS_FOLDER, str(self.wrID))

    def run(self):
        return NotImplementedError("Abstract class")

    def parse_call(self):
        module = importlib.import_module(self.context['package'])
        components = self.context['name'].split(".")
        for attr in components:
            module = getattr(module, attr)
        params = self.context['params'].copy()
        for key, value in self.context['params'].items():
            if value == "" or value == 0:
                del params[key]
        return module, params, components

    def log(self):
        """ Print the log of the finished process
        """
        param_str = []
        params = self.context['params']
        params.pop('filename', None)
        for k, v in params.items():
            v = "\'" + v + "\'" if type(v) == str else str(v)
            param_str.append(str(k) + "=" + v)
        params = ", ".join(param_str)

        call = self.context['package'] + "." + self.context['name']
        log = Process.objects.get(wrid=self.wrID, index=self.index)
        log.status = self.status
        log.output = self.output
        log.call = f"{call}(target, {params})"
        log.save()


class ReadStep(WorkerStep):

    def __init__(self, context, wrid, index, file, annData, subset=None):
        WorkerStep.__init__(self, context, wrid, index, file, annData)
        self.subset = subset

    def run(self):
        module, params, _ = self.parse_call()
        del params['filename']
        try:
            self.annData = module(self.file, **params)
            if self.subset is not None:
                self.annData = self.annData[self.subset, :]
        except Exception as e:
            self.output = str(e)
            self.status = 2
            self.log()
            return
        self.output = str(self.annData)
        self.status = 1
        self.log()


class ProcessStep(WorkerStep):
    def run(self):
        module, params, _ = self.parse_call()
        self.context['params'] = params
        try:
            module(self.annData, **params)
        except Exception:
            self.output = traceback.print_exc()
            self.status = 2
            self.log()
            return

        if self.context.get('view', None):
            self.annData.write(os.path.join(self.folder, f'views_{self.index}.h5ad'))

        self.output = str(self.annData)
        self.status = 1
        self.log()


class PlotStep(WorkerStep):
    def parse_call(self):
        module = importlib.import_module(self.context['package'])
        if self.context['package'] == "scanpy":
            module._settings.settings.figdir = self.folder
        components = self.context['name'].split(".")
        for attr in components:
            module = getattr(module, attr)
        params = self.context['params'].copy()
        for key, value in self.context['params'].items():
            if value == "" or value == 0:
                del params[key]
        return module, params, components

    def run(self):
        module, params, components = self.parse_call()
        self.context['params'] = params
        try:
            module(self.annData,
                   **params,
                   save=f'plot_{self.index}.png',
                   show=False)
        except Exception as e:
            self.output = str(e)
            self.status = 2
            self.log()
            return

        self.output = f"{components[-1]}plot_{self.index}.png"
        self.status = 1
        self.log()


class IPlotStep(WorkerStep):
    def run(self):
        module, params, components = self.parse_call()
        self.context['params'] = params
        try:
            with open(os.path.join(self.folder, f"{components[-1]}plot_{self.index}.json"), "w") as f:
                json.dump(module(self.annData,
                                 **params,
                                 save=os.path.join(self.folder, f"{components[-1]}plot_{self.index}.png")),
                          f)
        except Exception as e:
            self.output = str(e)
            self.status = 2
            self.log()
            return

        self.output = f"{components[-1]}plot_{self.index}.png"
        self.status = 1
        self.log()




class PredictStep(WorkerStep):

def run(self):
    module, params, _ = self.parse_call()
    self.context['params'] = params
    
    try:
        labels = model_run.run(query_dataloader)
        
    except Exception:
        self.output = traceback.print_exc()
        self.status = 2
        self.log()
        return

    if self.context.get('view', None):
        self.annData.write(os.path.join(self.folder, f'views_{self.index}.h5ad'))


    self.output = str(self.annData)
    self.status = 1
    self.log()
