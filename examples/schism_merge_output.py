#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Sample script to merge the output files from a number of simulations using the schism 
module. 

@author: khan
@email: jamal.khan@legos.obs-mip.fr
"""
import sys
sys.path.append('/home/khan/MEGA/Codes/SCHISMMB')

import glob
import os
import tarfile

from schism import io as schismio

inpath = '/run/media/khan/Storehouse/Projects/201803_Surge Level Under Current Climate/Experiments/Sensitivity/'
exps = glob.glob(os.path.join(inpath, 'Test*'))

for exp in exps:
    print(exp)
    path = os.path.join(exp, "outputs")
    os.mkdir(path)
    tar = tarfile.open(os.path.join(exp, "outputs_hydro.tar.gz"), mode='r:gz')
    tar.extractall(path)
    
    local2globals = schismio.Local2Globals(path)
    local2globals.load_files()
    local2globals.merge_nodes()

    nc = schismio.Schout(path=path, local2globals=local2globals, outpath=exp)
    nc.list_inputs()
    nc.create_file()
    for i in glob.glob(os.path.join(path, '*')):
        os.remove(i)
    os.rmdir(path)