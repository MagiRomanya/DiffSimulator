#!/usr/bin/env python3

import numpy as np
import csv


class SimulationRecorder:
    def __init__(self, filename="record.csv"):
        self.filename = filename
        self.wfile = open(filename, 'w')
        self.writer = csv.writer(self.wfile, delimiter=',',
                                 quoting=csv.QUOTE_NONNUMERIC)

    def __del__(self):
        self.wfile.close()

    def record_timestep(self, x: np.array, v: np.array):
        result = np.zeros(x.size*2)
        for i in range(x.size):
            result[i] = x[i]
            result[x.size+i] = v[i]
        self.writer.writerow(result)


class SimulationReader:
    def __init__(self, nDoF, filename="record.csv"):
        self.nDoF = nDoF
        self.filename = filename
        self.rfile = open(filename, 'r')
        self.reader = csv.reader(self.rfile, delimiter=',',
                                 quoting=csv.QUOTE_NONNUMERIC)

    def get_next_state(self):
        nDoF = self.nDoF
        row = self.reader.__next__()
        x = np.array(row[:nDoF])
        v = np.array(row[nDoF:])
        return (x, v)

    def from_the_start(self):
        del self.reader
        self.rfile.close()
        self.rfile = open(self.filename, 'r')
        self.reader = csv.reader(self.rfile, delimiter=',',
                                 quoting=csv.QUOTE_NONNUMERIC)

    def get_all_history(self):
        self.from_the_start()
        data = list(self.reader)
        return list(np.array(data).T)
