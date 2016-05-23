#!/usr/bin/env python

import os
import optparse
import re
import subprocess
from subprocess import *

import Bio
from Bio import *
from Bio.PDB import PDBParser, PDBIO

import rosetta
from rosetta import *
from rosetta.protocols.rigid import *
from rosetta.core.pack.rotamer_set import *
from rosetta.core.pack.interaction_graph import *
from rosetta.core.pack import *
from rosetta.core.scoring import *
from rosetta.core.graph import *
from toolbox import *
import StringIO

one_to_three = {'A': 'ALA',
                'R': 'ARG',
                'N': 'ASN',
                'D': 'ASP',
                'C': 'CYS',
                'E': 'GLU',
                'Q': 'GLN',
                'G': 'GLY',
                'H': 'HIS',
                'I': 'ILE',
                'L': 'LEU',
                'K': 'LYS',
                'M': 'MET',
                'F': 'PHE',
                'P': 'PRO',
                'S': 'SER',
                'T': 'THR',
                'W': 'TRP',
                'Y': 'TYR',
                'V': 'VAL',
            }
    
def rot_trans(pose, partners, flexibles, translation, rotation , jobs, out, mut):
 
    starting_p = pose
    dock_jump = 1

    setup_foldtree(starting_p, partners, Vector1([1]))

    scorefxn_talaris = create_score_function('talaris2014')
    
    dock_pert = RigidBodyPerturbMover(dock_jump, translation, rotation)
    
    movemap = MoveMap()
    movemap.set_jump(1, True)
    movemap.set_bb(True)
    
    tolerance = 0.01
    min_type = "dfpmin"
    minmover = MinMover(movemap, scorefxn_talaris, min_type, tolerance, True) 
    
    perturb = SequenceMover()
    perturb.add_mover(dock_pert)
    perturb.add_mover(minmover) 
    
    jd = PyJobDistributor(out+'/PDB/'+mut, jobs, scorefxn_talaris)
    jd.native_pose = starting_p
    counter=1
    while not jd.job_complete:
      starting_p=jd.native_pose
      perturb.apply(starting_p)
      
      compute_interactions(starting_p,'full.resfile', out+"/MAT/"+mut+"_"+str(counter)+".mat")
      command=["Convertor",out+"/MAT/"+mut+"_"+str(counter)+".mat",out+"/LG/"+mut+"_"+str(counter)+".LG"]
      call(command)
      
      command=["toulbar2",out+"/LG/"+mut+"_"+str(counter)+".LG"]
      tb2out=check_output(command)
      tb2out=tb2out.split('\n')
      for line in tb2out:
        line_split=line.split()
        if ("Optimum:" in line_split) and ("Energy:" in line_split):
            OptEnergy=float(line_split[3])
        elif ("Optimum:" in line_split) :
            OptSolution=line_split[1].split('-')
      OptSolution = [int(i) for i in OptSolution]
      
      get_Z_matrix(starting_p, OptSolution,OptEnergy, "full.resfile", flexibles,out+"/MAT/"+mut+"_"+str(counter)+".Zmat")
      command=["Convertor",out+"/MAT/"+mut+"_"+str(counter)+".Zmat",out+"/LG/"+mut+"_"+str(counter)+".ZLG"]
      call(command)
      
      jd.output_decoy(starting_p)
      counter += 1

def compute_interactions(pose,resfile,pdb_out):
    score_fxn = create_score_function('talaris2014')
    task_design = TaskFactory.create_packer_task(pose)
    task_design.initialize_from_command_line()
    parse_resfile(pose, task_design, resfile)
    pose.update_residue_neighbors()
    png = create_packer_graph(pose, score_fxn, task_design)  #Uncomment for latest Pyrosetta versions
    rotsets = RotamerSets()
    ig = pack_rotamers_setup(pose, score_fxn, task_design, rotsets)
    ig = InteractionGraphFactory.create_and_initialize_two_body_interaction_graph(task_design, rotsets, pose, score_fxn, png)  #Uncomment for latest Pyrosetta versions
    f = open(pdb_out,'w')
    ener=0.0
    f.write("Etemp = 0\n")
    for res1 in range(1,ig.get_num_nodes()+1):
        rotn=0
        prevresname=""
        for i in range(1, rotsets.rotamer_set_for_moltenresidue(res1).num_rotamers()+1):
            ener=str(ig.get_one_body_energy_for_node_state(res1,i))
            pres1 = rotsets.moltenres_2_resid(res1)
            resname = rotsets.rotamer_set_for_moltenresidue(res1).rotamer(i).name3()
            f.write("ROT::BKB\t"+str(res1)+"-"+resname+"-"+str(i)+" "+ener+" "+ener+'\n')
            prevresname=resname
    for res1 in range(1,ig.get_num_nodes()+1):
        for res2 in range(res1+1,ig.get_num_nodes()+1):
            if (ig.get_edge_exists(res1, res2)):
                rotn1=0
                prevres1name=""
                pres1 = rotsets.moltenres_2_resid(res1)
                for i in range(1, rotsets.rotamer_set_for_moltenresidue(res1).num_rotamers()+1):
                    nres1=rotsets.rotamer_set_for_moltenresidue(res1).rotamer(i).name3()
                    prevres1name=nres1
                    rotn2=0
                    prevres2name=""
                    pres2 = rotsets.moltenres_2_resid(res2)
                    for j in range(1, rotsets.rotamer_set_for_moltenresidue(res2).num_rotamers()+1):
                        ener=str(ig.get_two_body_energy_for_edge(res1,res2,i,j))
                        nres2=rotsets.rotamer_set_for_moltenresidue(res2).rotamer(j).name3()
                        f.write("ROT::ROT\t"+str(res1)+"-"+nres1+"-"+str(i)+"::"+str(res2)+"-"+nres2+"-"+str(j)+" "+ener+" "+ener+'\n')
                        prevres2name=nres2

    f.close()
    return pose
    
def get_Z_matrix(pose, optsolution, optenergy, resfile, flexibles,out_matrix):
    score_fxn = create_score_function('talaris2014')
    task_design = TaskFactory.create_packer_task(pose)
    task_design.initialize_from_command_line()
    parse_resfile(pose, task_design, resfile)
    pose.update_residue_neighbors()
    png = create_packer_graph(pose, score_fxn, task_design)  #Uncomment for latest Pyrosetta versions
    rotsets = RotamerSets()
    ig = pack_rotamers_setup(pose, score_fxn, task_design, rotsets)
    ig = InteractionGraphFactory.create_and_initialize_two_body_interaction_graph(task_design, rotsets, pose, score_fxn, png)  #Uncomment for latest Pyrosetta versions
    mat=optsolution
    template_energy = optenergy
    for i in range(0, len(flexibles)):
        template_energy-=ig.get_one_body_energy_for_node_state(flexibles[i],int(mat[flexibles[i]-1]+1))
#        print str(flexibles[i])+" "+str(ig.get_one_body_energy_for_node_state(flexibles[i],int(mat[flexibles[i]-1]+1)))
#        print "Removing "+str(ig.get_one_body_energy_for_node_state(flexibles[i],int(mat[flexibles[i]-1]+1)))+" from unary term "+str(flexibles[i])
        for j in range(0,len(mat)):
            if (ig.get_edge_exists(flexibles[i], j+1)):
                if (j+1) in flexibles:
                    if (flexibles[i]<j+1):
                        template_energy-=ig.get_two_body_energy_for_edge(flexibles[i],j+1,int(mat[flexibles[i]-1]+1),int(mat[j]+1))
#                        print str(flexibles[i])+" "+str(j+1)+" "+str(ig.get_two_body_energy_for_edge(flexibles[i],j+1,int(mat[flexibles[i]-1]+1),int(mat[j]+1)))+" "+str(ig.get_two_body_energy_for_edge(j+1,flexibles[i],int(mat[j]+1),int(mat[flexibles[i]-1]+1)))
                else:
                    if (flexibles[i]<j+1):
                        template_energy-=ig.get_two_body_energy_for_edge(flexibles[i],j+1,int(mat[flexibles[i]-1]+1),int(mat[j]+1))
                    else:
                        template_energy-=ig.get_two_body_energy_for_edge(j+1,flexibles[i],int(mat[j]+1),int(mat[flexibles[i]-1]+1))
 #                   print str(flexibles[i])+" "+str(j+1)+" "+str(ig.get_two_body_energy_for_edge(flexibles[i],j+1,int(mat[flexibles[i]-1]+1),int(mat[j]+1)))+" "+str(ig.get_two_body_energy_for_edge(j+1,flexibles[i],int(mat[j]+1),int(mat[flexibles[i]-1]+1)))
 #               print "Removing "+str(ig.get_two_body_energy_for_edge(flexibles[i],j+1,int(mat[flexibles[i]-1]+1),int(mat[j]+1)))+" from binary term "+str(flexibles[i])+" "+str(j+1)
#                ig.clear_two_body_energy_for_edge(flexibles[i],j+1,int(mat[flexibles[i]-1]+1),int(mat[j]+1))
    
    f = open(out_matrix,'w')
    f.write("Etemp = "+str(template_energy)+"\n")
    binary_terms=StringIO.StringIO()
    for res1 in range(0, len(flexibles)):
        for i in range(1, rotsets.rotamer_set_for_moltenresidue(flexibles[res1]).num_rotamers()+1):
            nres1=rotsets.rotamer_set_for_moltenresidue(flexibles[res1]).rotamer(i).name3()
            unary_ener=ig.get_one_body_energy_for_node_state(flexibles[res1],i)
            for res2 in range(1,ig.get_num_nodes()+1):
                if (ig.get_edge_exists(flexibles[res1], res2)):
                    if res2 in flexibles:
                        if (flexibles[res1]<res2):
                            for j in range(1, rotsets.rotamer_set_for_moltenresidue(res2).num_rotamers()+1):
                                ener=str(ig.get_two_body_energy_for_edge(flexibles[res1],res2,i,j))
                                nres2=rotsets.rotamer_set_for_moltenresidue(res2).rotamer(j).name3()
                                binary_terms.write("ROT::ROT\t"+str(flexibles[res1])+"-"+nres1+"-"+str(i)+"::"+str(res2)+"-"+nres2+"-"+str(j)+" "+ener+" "+ener+'\n')
                    else:
                        unary_ener+=ig.get_two_body_energy_for_edge(flexibles[res1],res2,i,mat[res2-1]+1)
            f.write("ROT::BKB\t"+str(flexibles[res1])+"-"+nres1+"-"+str(i)+" "+str(unary_ener)+" "+str(unary_ener)+'\n')
    f.write(binary_terms.getvalue())
    binary_terms.close()
    f.close()

def mutation_rot_trans(pdb_file, seq_file, jobs):
  if os.path.exists( os.getcwd() + '/' + pdb_file ) and pdb_file:
    init()
    pose=Pose()
    pose=pose_from_pdb(pdb_file)
    input_file_name=os.getcwd() + '/' + pdb_file.split('.pdb')[0]

    ## create a score function
    scorefxn = create_score_function('talaris2014')
    
    ## parse the fasta sequence file into a dictionary : sequences[#position]=[list of mutation]
    sequences={}
    sequences[pdb_file.split('.pdb')[0]]=[]
    if os.path.exists( os.getcwd() + '/' + seq_file ) and seq_file:
        r = re.compile("([0-9]+)([a-zA-Z]+)")
        seqfile=open(seq_file,'r')
        for seq in seqfile.readlines():
            sequences[seq[:-1]]=[r.match(i).groups() for i in seq.split('_')]
            
    ## parse the flexible file to get the flexibles residues
    flexibles_rec=open("flexibles.receptor",'r')
    flexibles_rec=flexibles_rec.readlines()[0]
    flexibles_rec=[int(i) for i in flexibles_rec.split()]
    
    flexibles_lig=open("flexibles.ligand",'r')
    flexibles_lig=flexibles_lig.readlines()[0]
    flexibles_lig=[int(i) for i in flexibles_lig.split()]
    
    flexibles=sorted(flexibles_rec+flexibles_lig)
    
    ## First minimisation (may do fastrelax ?) 
    movemap = MoveMap()
    movemap.set_jump(1, True)
    movemap.set_bb(True)
    tolerance = 0.01
    min_type = "dfpmin"
    minmover = MinMover(movemap, scorefxn, min_type, tolerance, True) 
    minmover.apply(pose)
    
    pose.dump_pdb(input_file_name+"_min.pdb")
    
    ###### Mutation Loop ####
    for mut in sequences.keys():
      mut_pose=pose
      print "Processing Mutation:",mut
      
      for mut_tuple in sequences[mut]:
        index = int(mut_tuple[0]) + 1
        aa = mut_tuple[1]
        mutator=MutateResidue(int(index), one_to_three[aa])
        mutator.apply(mut_pose)
        
      mut_folder=os.getcwd()+"/"+mut
      if not os.path.exists(mut_folder):
        os.mkdir(mut_folder)
        
      if not os.path.exists(mut_folder+'/LG'):
        os.mkdir(mut_folder+'/LG')

      if not os.path.exists(mut_folder+'/PDB'):
        os.mkdir(mut_folder+'/PDB')
        
      if not os.path.exists(mut_folder+'/MAT'):
        os.mkdir(mut_folder+'/MAT')
        
      mut_pose.dump_pdb(mut_folder+"/PDB/"+mut+'.pdb')

      #~ pose_packer = standard_packer_task(mut_pose)
      #~ pose_packer.restrict_to_repacking()
      #~ packmover = PackRotamersMover(scorefxn, pose_packer)
      #~ packmover.apply(mut_pose)
        
      ## Minimization on the mutant pose. Useless ?
      minmover.apply(mut_pose)
      mut_pose.dump_pdb(mut_folder+"/PDB/"+mut+'_min.pdb')

      io = PDBIO()
      pdb = PDBParser().get_structure(mut, mut_folder+"/PDB/"+mut+"_min.pdb")
      chain_name=[]
      chain_length=[]
      for chain in pdb.get_chains():
        io.set_structure(chain)
        io.save(mut_folder+'/PDB/'+pdb.get_id() + "_" + chain.get_id() + ".pdb")
        chain_length.append(len(chain))
        chain_name.append(chain.get_id())
      
      ## Separate the two chains ex: E_I or A_B. Need to rename if more than two chains ?
      pose_prot_1=Pose()
      pose_prot_2=Pose()
      pose_prot_1=pose_from_pdb(mut_folder+'/PDB/'+pdb.get_id() + "_" + chain_name[0] + ".pdb")
      pose_prot_2=pose_from_pdb(mut_folder+'/PDB/'+pdb.get_id() + "_" + chain_name[1] + ".pdb")
      #repack a faire ???? 
      
      ## Minimise the lonely partners
      minmover.apply(pose_prot_1)
      minmover.apply(pose_prot_2)
      
      ###### Compute FULL SCP matrix, Calculate the optimal solution and optimal energy and compute Z matrix.
      ##### FOR THE RECEPTOR
      compute_interactions(pose_prot_1,'full.resfile', mut_folder+'/MAT/receptor.mat')
      command=["Convertor",mut_folder+'/MAT/receptor.mat',mut_folder+'/LG/receptor.LG']
      call(command)
      
      command=["toulbar2",mut_folder+"/LG/receptor.LG"]
      tb2out=check_output(command)
      tb2out=tb2out.split('\n')
      OptEnergy=0
      
      for line in tb2out:
        line_split=line.split()
        if ("Optimum:" in line_split) and ("Energy:" in line_split):
            OptEnergy=float(line_split[3])
        elif ("Optimum:" in line_split) :
            OptSolution=line_split[1].split('-')
      OptSolution = [int(i) for i in OptSolution]
      
      get_Z_matrix(pose_prot_1, OptSolution, OptEnergy, "full.resfile", flexibles_rec, mut_folder+"/MAT/receptor.Zmat")	
      command=["Convertor",mut_folder+"/MAT/receptor.Zmat",mut_folder+"/LG/receptor.ZLG"]
      call(command)
      
      ##### FOR THE LIGAND
      flexibles_lig_renum=[i-int(chain_length[0]) for i in flexibles_lig]
      compute_interactions(pose_prot_2,'full.resfile', mut_folder+'/MAT/ligand.mat')
      command=["Convertor",mut_folder+"/MAT/ligand.mat",mut_folder+"/LG/ligand.LG"]
      call(command)
      
      command=["toulbar2",mut_folder+"/LG/ligand.LG"]
      tb2out=check_output(command)
      tb2out=tb2out.split('\n')
      for line in tb2out:
        line_split=line.split()
        if ("Optimum:" in line_split) and ("Energy:" in line_split):
            OptEnergy=float(line_split[3])
        elif ("Optimum:" in line_split) :
            OptSolution=line_split[1].split('-')
      OptSolution = [int(i) for i in OptSolution]
      
      get_Z_matrix(pose_prot_2,OptSolution,OptEnergy,"full.resfile",flexibles_lig_renum,mut_folder+"/MAT/ligand.Zmat")		
      command=["Convertor",mut_folder+"/MAT/ligand.Zmat",mut_folder+"/LG/ligand.ZLG"]
      call(command)
     
      # FOR THE COMPLEX
      compute_interactions(mut_pose,'full.resfile', mut_folder+"/MAT/"+mut+'_min.mat')
      command=["Convertor",mut_folder+"/MAT/"+mut+'_min.mat',mut_folder+"/LG/"+mut+'_min.LG']
      call(command)
      
      command=["toulbar2",mut_folder+"/LG/"+mut+'_min.LG']
      tb2out=check_output(command)
      tb2out=tb2out.split('\n')
      for line in tb2out:
        line_split=line.split()
        if ("Optimum:" in line_split) and ("Energy:" in line_split):
            OptEnergy=float(line_split[3])
        elif ("Optimum:" in line_split) :
            OptSolution=line_split[1].split('-')
      OptSolution = [int(i) for i in OptSolution]
      
      get_Z_matrix(mut_pose,OptSolution,OptEnergy,"full.resfile",flexibles,mut_folder+"/MAT/"+mut+'_min.Zmat')
      command=["Convertor",mut_folder+"/MAT/"+mut+'_min.Zmat',mut_folder+"/LG/"+mut+'_min.ZLG']
      call(command)
      
      partners=chain_name[0]+'_'+chain_name[1]
      rot_trans(mut_pose, partners, flexibles, 1, 3, jobs, mut_folder,mut)

      print "Finish Processing Mutation:",mut
  else:
    "ERROR: PDB FILE NOT EXISTING"

        

parser=optparse.OptionParser()
parser.add_option('--pdb', dest = 'pdb_file',
    default = '',    
    help = 'the backbone in PDB format' )

parser.add_option('--seq', dest = 'seq_file',
    default = '',    
    help = 'the sequences to map' )

parser.add_option( '--jobs', dest='jobs' ,
    default = '1',    # default to single trajectory for speed
    help = 'the number of jobs (trajectories) to perform')
    
    
(options,args) = parser.parse_args()

(options,args) = parser.parse_args()

pdb_file=options.pdb_file
sequence_file = options.seq_file
jobs = int(options.jobs)


################# MUTATION, PDB and MATRIX PRODUCTION #############
        
mutation_rot_trans(pdb_file,sequence_file, jobs)

