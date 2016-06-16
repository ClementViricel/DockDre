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

parser.add_option('--nb', dest = 'nb',
    default = 5,    
    help = 'Number of docked pose taken as instance in Z computation' )

(options,args) = parser.parse_args()

pdb_file=options.pdb_file
seq_file = options.seq_file
nb_instance=options.nb

sequences=[]
sequences.append(pdb_file.split('.pdb')[0])

## parse the flexible file to get the flexibles residues
flexibles_rec=open("flexibles.receptor",'r')
flexibles_rec=flexibles_rec.readlines()[0]
flexibles_rec=[int(i) for i in flexibles_rec.split()]

flexibles_lig=open("flexibles.ligand",'r')
flexibles_lig=flexibles_lig.readlines()[0]
flexibles_lig=[int(i) for i in flexibles_lig.split()]

flexibles=sorted(flexibles_rec+flexibles_lig)


if os.path.exists( os.getcwd() + '/' + seq_file ) and seq_file:
  seqfile=open(seq_file,'r')
  for seq in seqfile.readlines():
    sequences.append(seq[:-1])
  seqfile.close()
  
  K=[]
  time=[]
  for mut in sequences:
    print mut
    if os.path.exists( os.getcwd() + '/' + mut ):
      Kmut=0
      LogZb=[]
      out=open(os.getcwd()+'/'+mut+'/K_list','w')
      scorefile=open(os.getcwd()+'/'+mut+'/score.sc','r')
      score_file=scorefile.readlines()
      for instance in score_file[0:nb_instance]:
        instance=instance.split()
        print "Pose:", instance[0]
        out.write(str(instance[0])+' ')
        f=open(os.getcwd()+'/'+mut+"/SOL/"+mut+"_"+instance[0]+".sol",'r')
        sol=[int(i) for i in f.readlines()[0].split()]
        f.close()
        opt=[]
        for i in range(0,len(sol)):
          if not i+1 in flexibles:
            opt+=","+str(i)+"="+str(sol[i])
        command=["toulbar2",os.getcwd()+'/'+mut+"/LG/"+mut+"_"+instance[0]+".LG","-logz","-ztmp","-x="+opt]
        tb2out=check_output(command)
        tb2out=tb2out.split('\n')
        for line in tb2out:
          if "<= Log(Z) <=" in line:
           LogZb.append(float(line.split()[0]))
           out.write(line.split()[0]+' '+line.split()[12]+'\n')
           time.append(float(line.split()[12]))
      LogZ=logsumexp(LogZb)
      
      print "Receptor"
      f=open(os.getcwd()+'/'+mut+"/SOL/receptor.sol",'r')
      sol_rec=[int(i) for i in f.readlines()[0].split()]
      f.close()
      opt=[]
      for i in range(0,len(sol_rec)):
        if not i+1 in flexibles_rec:
          opt+=","+str(i)+"="+str(sol_rec[i])
      command=["toulbar2",os.getcwd()+'/'+mut+"/LG/receptor.LG","-logz","-ztmp","-x="+opt]
      tb2out=check_output(command)
      tb2out=tb2out.split('\n')
      for line in tb2out:
        if "<= Log(Z) <=" in line:
          LogZrec=float(line.split()[0])
          time.append(float(line.split()[12]))
          out.write("Receptor "+line.split()[0]+' '+line.split()[12]+'\n')
          
      print "Ligand"
      flexibles_lig_renum=[i-int(len(sol_rec)) for i in flexibles_lig]
      f=open(os.getcwd()+'/'+mut+"/SOL/ligand.sol",'r')
      sol_lig=[int(i) for i in f.readlines()[0].split()]
      f.close()
      opt=[]
      for i in range(0,len(sol_lig)):
        if not i+1 in flexibles_lig_renum:
          opt+=","+str(i)+"="+str(sol_lig[i])
      command=["toulbar2",os.getcwd()+'/'+mut+"/LG/ligand.LG","-logz","-ztmp","-x="+opt]
      tb2out=check_output(command)
      tb2out=tb2out.split('\n')
      for line in tb2out:
        if "<= Log(Z) <=" in line:
          LogZlig=float(line.split()[0])
          time.append(float(line.split()[12]))
          out.write("Ligand "+line.split()[0]+' '+line.split()[12]+'\n')
          
      Kmut=LogZ-LogZrec-LogZlig
      scorefile.close()
      out.close()
      K.append([mut,Kmut])
  print K
  K=sorted(K,key=itemgetter(1))
  time=sum(time)
  out_file=open(os.getcwd()+"/K_list",'w')
  for k in K:
    out_file.write(k[0]+' '+str(k[1])+'\n')
  out_file.close()
  time_file=open(os.getcwd()+"/time.log",'w')
  time_file.write(str(time))
  time_file.close()
    

    
