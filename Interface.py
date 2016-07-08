#!/usr/bin/env python
from datetime import datetime
import os
import optparse
import re
import subprocess
from subprocess import *
import numpy as np
from operator import sub

import Bio
from Bio import *
from Bio.PDB import *

import rosetta
from rosetta import *
from rosetta.protocols.rigid import *
from rosetta.core.pack.rotamer_set import *
from rosetta.core.pack.interaction_graph import *
from rosetta.core.pack import *
from rosetta.core.scoring import *
from rosetta.core.graph import *
from rosetta.protocols.scoring import Interface 

from toolbox import *
import StringIO

parser=optparse.OptionParser()
parser.add_option('--pdb', dest = 'pdb_file', default = '', help = 'Protein comple in PDB format' )
parser.add_option( '--dist', dest='dist' ,
    default = 10.0,   
    help = 'Size of interface')

(options,args) = parser.parse_args()
pdb_file=options.pdb_file
dist=float(options.dist)

init()
pose=pose_from_pdb(pdb_file)
setup_foldtree(pose, "A_B", Vector1([1]))
score=create_score_function("talaris2014")
score(pose)
interface = Interface(1)
interface.distance(dist)
interface.calculate(pose)
contact = interface.contact_list()
interface_residue = interface.pair_list()

print pose.fold_tree()
print interface_residue
