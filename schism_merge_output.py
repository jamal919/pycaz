# -*- coding: utf-8 -*-
"""
Created on Tue Jul 24 14:56:03 2018

@author: khan
"""
import glob
import os
import tarfile
from SCHISM import schismio

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