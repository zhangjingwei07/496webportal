import json
import os
from threading import Thread

from dataset.models import DataSet
from dataset.utils import get_anndata_attrs
from settings.settings import USER_PROCESS_FOLDER
from .models import Process, WorkerRecord
from .worker_step import ReadStep, ProcessStep, PlotStep, IPlotStep


class Worker(Thread):
    """ A process worker that runs and stores the resulted of the defined procedure
        @:param process: The list of dictionary contains the whole process
        @:param curr: the index of the current running function
        @:param status: -1 suspended, 0 pending, 1 running, 2 finished
        @:param id: the unique id indicating a running worker, used to retrieve
                   the worker
        @:param name: the name of this worker

    """

    def __init__(self, process, name):
        Thread.__init__(self)
        self.process = process
        self.curr = 0
        self.status = 0
        self.id = ""
        self.annData = None
        self.filename = None
        self.name = name

    def check_integrity(self):
        """ Check whether this process can run
        """

        # file-check
        reader = self.process[0]
        file_id = reader['params'].get('filename', None)
        if not file_id:
            return "Missing File", False
        self.filename = DataSet.objects.get(id=file_id).path
        if not self.filename:
            return {'info': "Missing File", 'status': False}

        wr = WorkerRecord(status=0,
                          curr=0,
                          total=len(self.process),
                          name=self.name)
        wr.save()
        self.id = wr.id

        path = os.path.join(USER_PROCESS_FOLDER, str(self.id))
        print(path)
        os.mkdir(path)

        with open(os.path.join(path, "procedure.json"), "w") as f:
            json.dump(self.process, f)

        index = 0
        for p in self.process:
            process = Process(wrid=self.id,
                              index=index,
                              call=p['package'] + "." + p['name'],
                              status=0,
                              output="",
                              type=p['type'])
            index += 1

            process.save()

        return {'info': self.name, 'status': True}

    def run(self):
        wr = WorkerRecord.objects.get(id=self.id)
        for step in self.process:
            if step['type'] == 'reader':
                worker_step = ReadStep(step, self.id, self.curr, self.filename, None)
            elif step['type'] == 'processing':
                worker_step = ProcessStep(step, self.id, self.curr, "", self.annData)
            elif step['type'] == 'plot':
                worker_step = PlotStep(step, self.id, self.curr, "", self.annData)
            elif step['type'] == 'iplot':
                worker_step = IPlotStep(step, self.id, self.curr, "", self.annData)
            else:
                wr.status = 2
                wr.save()
                return
            worker_step.run()

            if worker_step.status == 2:
                wr.status = 2
                wr.save()
                return
            elif worker_step.status == 1:
                self.annData = worker_step.annData
                self.curr += 1
                wr.curr = self.curr
                wr.save()
        wr.status = 1
        wr.save()
        path = os.path.join(USER_PROCESS_FOLDER, str(self.id), "results.h5ad")
        self.annData.write(path)

        saved_file = DataSet(
            user='__RESULT',
            name=f'Worker_{self.id}',
            path=path,
            description="",
            n_obs=self.annData.n_obs,
            n_vars=self.annData.n_vars,
            attrs=json.dumps(get_anndata_attrs(self.annData))
        )
        saved_file.save()
        return self.id


def log_sync(adata, wrid, call, type="processing"):
    wr = WorkerRecord.objects.get(id=wrid)
    index = wr.total + 1
    wr.total = index
    wr.curr = index
    wr.save()
    process = Process(wrid=wrid,
                      index=index,
                      call=call,
                      status=1,
                      output=str(adata),
                      type=type)
    process.save()
    saved_file = DataSet.objects.get(user="__RESULT", name=f'Worker_{wrid}')
    saved_file.n_obs = adata.n_obs
    saved_file.n_vars = adata.n_vars
    saved_file.attrs = json.dumps(get_anndata_attrs(adata))
    saved_file.save()
    adata.write(saved_file.path)
