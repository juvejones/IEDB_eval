#!/usr/bin/python

##Some universal classes across the module

class PepData(object):
    def __init__(self,size):
        self.name = [None]*size
        self.meas_nm = [None]*size
        self.meas_bi = [None]*size
        self.meas_contin = [None]*size
        self.predict = [None]*size
        self.predict_rank = [None]*size

    def __enter__(self):
        return self

    def close(self):
    	print("")

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()


class Dir(object):
	def __init__(self):
		self.logfile = ""
		self.input_dir = ""
		self.tmp_dir = ""
		self.output_dir = ""
		self.data_dir = ""