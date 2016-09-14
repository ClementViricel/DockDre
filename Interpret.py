#!/usr/bin/env python

import os
import optparse
import re
import subprocess
from subprocess import *
import numpy as np
from operator import *
from scipy.misc import logsumexp


parser=optparse.OptionParser()
parser.add_option('--pdb', dest = 'pdb_file',
    default = '',    
    help = 'Protein comple in PDB format' )
    
parser.add_option('--seq', dest = 'seq_file',
    default = '',    
    help = 'Sequences file' )

(options,args) = parser.parse_args()

pdb_file=options.pdb_file
seq_file = options.seq_file

sequences=[]
sequences.append(pdb_file.split('.pdb')[0])


if os.path.exists( os.getcwd() + '/' + seq_file ) and seq_file:
  seqfile=open(seq_file,'r')
  for seq in seqfile.readlines():
    sequences.append(seq[:-1])
  seqfile.close()
  K=[]
  time=[]
  for mut in sequences:
    if os.path.exists( os.getcwd() + '/' + mut ):
      Kmut=0
      LogZb=[]
      files=os.listdir(os.getcwd() + '/' + mut + '/Z')
      out=open(os.getcwd()+'/'+mut+'/K_list','w')
      for file_Z in files:
        if file_Z=="receptor.Zlog":
          f=open(os.getcwd() + '/' + mut + '/Z/'+file_Z,'r')
          for line in f:
            if "<= Log(Z) <=" in line:
              LogZrec=float(line.split()[0])
              time.append(float(line.split()[12]))
              out.write("Receptor "+line.split()[0]+' '+line.split()[12]+'\n')

        elif file_Z=="ligand.Zlog":
          f=open(os.getcwd() + '/' + mut + '/Z/'+file_Z,'r')
          for line in f:
            if "<= Log(Z) <=" in line:
              LogZlig=float(line.split()[0])
              time.append(float(line.split()[12]))
              out.write("Ligand "+line.split()[0]+' '+line.split()[12]+'\n')
        else:
          f=open(os.getcwd() + '/' + mut + '/Z/'+file_Z,'r')
          for line in f:
            if "<= Log(Z) <=" in line:
              LogZb.append(float(line.split()[0]))
              out.write(file_Z+" "+line.split()[0]+' '+line.split()[12]+'\n')
              time.append(float(line.split()[12]))
      LogZ=logsumexp(LogZb)
      Kmut=LogZ-LogZrec-LogZlig
      out.close()
      K.append([mut,Kmut])
  K=sorted(K,key=itemgetter(1))
  time=sum(time)
  out_file=open(os.getcwd()+"/K_list",'w')
  for k in K:
    out_file.write(k[0]+' '+str(k[1])+'\n')
  out_file.close()
  time_file=open(os.getcwd()+"/time.log",'w')
  time_file.write(str(time))
  time_file.close()
    

    
