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

reference_energy={
    'ALA': 0.773742,
    'ARG': 0.443793,
    'ASN': -1.63002,
    'ASP': -1.96094,
    'CYS': 0.61937,
    'CYD': 0.61937,
    'GLU':0.173326,
    'GLN': 0.388298,
    'GLY': 1.0806,
    'HIS': -0.358574,
    'ILE': 0.761128,
    'LEU': 0.249477,
    'LYS': -1.19118,
    'MET': -0.250485,
    'PHE': -1.51717,
    'PRO': -0.32436,
    'SER': 0.165383,
    'THR': 0.20134,
    'TRP': 0.979644,
    'TYR': 1.23413,
    'VAL': 0.162496,
    }



def compute_interactions(pose, resfile, out, score_fxn):

    ## Minimization
    movemap = MoveMap()
    movemap.set_jump(1, True)
    movemap.set_bb(True)
    tolerance = 0.001
    min_type = "dfpmin"
    minmover = MinMover(movemap, score_fxn, min_type, tolerance, True) 
    minmover.apply(pose)
    #relax=FastRelax(score_fxn)
    #relax.apply(pose)
    
    task_design = TaskFactory.create_packer_task(pose)
    task_design.initialize_from_command_line()
    parse_resfile(pose, task_design, resfile)
    pose.update_residue_neighbors()
    png = create_packer_graph(pose, score_fxn, task_design)  #Uncomment for latest Pyrosetta versions
    rotsets = RotamerSets()
    ig = pack_rotamers_setup(pose, score_fxn, task_design, rotsets)
    ig = InteractionGraphFactory.create_and_initialize_two_body_interaction_graph(task_design, rotsets, pose, score_fxn, png)  #Uncomment for latest Pyrosetta versions
    out_uai=out.split('.LG')[0]+'.uai'
    g = open(out,'w')
    h = open(out_uai,'w')
    ener=0.0
    g.write("MARKOV\n")
    g.write(str(ig.get_num_nodes())+'\n')
    h.write("MARKOV\n")
    h.write(str(ig.get_num_nodes())+'\n')
    domain=StringIO.StringIO()
    scope=StringIO.StringIO()
    unary_terms=StringIO.StringIO()
    exp_unary_terms=StringIO.StringIO()
    binary_terms=StringIO.StringIO()
    exp_binary_terms=StringIO.StringIO()
    number_of_functions=0
    for res1 in range(1,ig.get_num_nodes()+1):
        number_of_functions += 1
        domain_res=str(rotsets.rotamer_set_for_moltenresidue(res1).num_rotamers())
        domain.write(domain_res+' ')
        scope.write("1 "+str(res1-1)+'\n')
        unary_terms.write(domain_res+'\n')
        exp_unary_terms.write(domain_res+'\n')
        for i in range(1, rotsets.rotamer_set_for_moltenresidue(res1).num_rotamers()+1):
            resname = rotsets.rotamer_set_for_moltenresidue(res1).rotamer(i).name3()
            ener=ig.get_one_body_energy_for_node_state(res1,i) + reference_energy[resname]
            unary_terms.write(str(-ener)+' ')
            exp_unary_terms.write(str(np.exp(-ener))+' ')
        unary_terms.write('\n')
        exp_unary_terms.write('\n')
    domain.write('\n')

    for res1 in range(1,ig.get_num_nodes()+1):
        for res2 in range(res1+1,ig.get_num_nodes()+1):
            if (ig.get_edge_exists(res1, res2)):
                number_of_functions += 1
                scope.write("2 "+str(res1-1)+" "+str(res2-1)+"\n")
                domain_res=str(rotsets.rotamer_set_for_moltenresidue(res1).num_rotamers()*rotsets.rotamer_set_for_moltenresidue(res2).num_rotamers())
                binary_terms.write(domain_res+'\n')
                exp_binary_terms.write(domain_res+'\n')
                for i in range(1, rotsets.rotamer_set_for_moltenresidue(res1).num_rotamers()+1):
                    for j in range(1, rotsets.rotamer_set_for_moltenresidue(res2).num_rotamers()+1):
                        ener=ig.get_two_body_energy_for_edge(res1,res2,i,j)
                        binary_terms.write(str(-ener)+' ')
                        exp_binary_terms.write(str(np.exp(-ener))+' ')
                    binary_terms.write('\n')
                    exp_binary_terms.write('\n')

    g.write(domain.getvalue())
    h.write(domain.getvalue())
    g.write(str(number_of_functions)+'\n')
    h.write(str(number_of_functions)+'\n')
    g.write(scope.getvalue())
    h.write(scope.getvalue())
    
    g.write(unary_terms.getvalue())
    g.write(binary_terms.getvalue())

    h.write(exp_unary_terms.getvalue())
    h.write(exp_binary_terms.getvalue())

    domain.close()
    scope.close()
    unary_terms.close()
    binary_terms.close()
    exp_unary_terms.close()
    exp_binary_terms.close()
    g.close()
    h.close()



parser=optparse.OptionParser()

parser.add_option('--mut', dest = 'mut',
    default = '',    
    help = 'Sequences to map' )

parser.add_option( '--delta', dest='delta' ,
    default = 1.0,    
    help = 'Size of translation segment')
  
parser.add_option( '--teta', dest='teta' ,
    default = 3.0,   
    help = 'Size of rotation segment')

parser.add_option( '--out', dest='out' ,
    default = "out",   
    help = 'out file')

parser.add_option( '--count', dest='counter' ,
    default = "out",   
    help = 'out file')

parser.add_option( '--beta', dest='beta' ,
    default = False,   
    help = 'Beta scoring function')    

(options,args) = parser.parse_args()

delta=float(options.delta)
teta=float(options.teta)
mut=options.mut
out=options.out
counter=int(options.counter)
beta=options.beta

if beta == True:
    init('-ex1 false -ex1aro false -ex2 false -ex2aro false -beta')
    score_fxn = create_score_function('beta_july15')
elif beta == False:
    init('-ex1 false -ex1aro false -ex2 false -ex2aro false')
	score_fxn = create_score_function('talaris2014')
        
pose=Pose()
pose=pose_from_pdb(out+'/PDB/'+mut+'_'+str(counter)+"_("+str(delta)+")T_("+str(teta)+")R.pdb")
setup_foldtree(pose,"A_B",Vector1([1]))

compute_interactions(pose,'full.resfile', out+"/LG/"+mut+"_"+str(counter)+".LG",score_fxn) ## LG for FSCP
pose.dump_pdb(out+"/PDB/"+mut+"_"+str(counter)+"_min.pdb")
os.rename(out+"/LG/"+mut+"_"+str(counter)+".uai", out+"/UAI/"+mut+"_"+str(counter)+".uai")

command=["toulbar2",out+"/LG/"+mut+"_"+str(counter)+".LG","-w="+out+"/SOL/"+mut+"_"+str(counter)+".sol"] # FSCP
tb2out=check_output(command)
tb2out=tb2out.split('\n')

for line in tb2out:
  line_split=line.split()
  if ("Optimum:" in line_split) and ("Energy:" in line_split):
    OptEnergy=float(line_split[3])

f=open(out+"/score.sc",'ab')
f.write(str(counter)+' '+str(OptEnergy)+"\n")
f.close()
