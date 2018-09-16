# -*- coding: utf-8 -*-
"""
Reads tidefac output and adapt given bctides file as required by the tidefac 
outputs.

@author: khan
@email: jamal.khan@legos.obs-mip.fr
"""
from __future__ import print_function
import numpy as np
import re

class Bctides(object):
    def __init__(self):
        pass

class Tidefacout(object):
    def __init__(self, filepath):
        self.filepath = filepath
        self.r = re.compile(r'[0-9]?.[0-9]\d+')
        self.__read()
    
    def __read(self):
        # Reading date information
        with open(self.filepath, 'r') as f:
            # Reading the date section
            __ds = f.readline()
            __date = np.array(self.r.findall(__ds), dtype='float')
            self.hour = __date[0]
            self.day = __date[1]
            self.month = __date[2]
            self.year = __date[3]
            
            # Reading the run length section
            __ds = f.readline()
            __rnday = self.r.findall(__ds)
            
            if(len(__rnday)==0):
                __ds = f.readline()
                
            __rnday = np.array(self.r.findall(__ds), dtype='float')
            self.rnday = __rnday[0]

        # Reading the constants, node factor and eq. argument ref. to GM in deg.
        self.const = np.genfromtxt(fname=self.filepath, dtype=None, skiprows=8, \
                                delimiter=None, names=True, autostrip=True)
        self.const = np.array([[i for i in j] for j in self.const])
        print(self.const.shape)

    def __str__(self):
        __msg = 'Tidefac output for {:.1f} days run starting from {:4.0f}/{:02.0f}/{:02.0f} {:02.1f} UTC'.format(self.rnday,\
                self.year, self.month, self.day, self.hour)
        return(__msg)

if __name__=='__main__':
    file = '/run/media/khan/Workbench/Projects/Surge Model/Benchmark_Model_Mesh/Test02/tide_fac.out'
    tfout = Tidefacout(filepath=file)
    print(tfout)
