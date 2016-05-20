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
from rosetta.core.pack.rotamer_set import *
from rosetta.core.pack.interaction_graph import *
from rosetta.core.pack import *
from rosetta.core.scoring import *
from rosetta.core.graph import *
from toolbox import *

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
    
def rot_trans(pose, partners,
        translation = 0.2, rotation = 3.0,
        jobs = 1, out = 'dock_output'):
 
    starting_p = pose
    dock_jump = 1

    setup_foldtree(starting_p, partners, Vector1([1]))

    scorefxn_docking = create_score_function_ws_patch('talaris2014','docking')
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
    
    jd = PyJobDistributor(out, jobs, scorefxn_talaris)
    jd.native_pose = starting_p
    while not jd.job_complete:
      starting_p=jd.native_pose
      perturb.apply(starting_p)
      jd.output_decoy(starting_p)

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
    
def get_template_energy(pose, optsolution, optenergy, resfile, flexibles):
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
    print len(mat)
    template_energy = optenergy
    for i in range(0, len(flexibles)):
        print int(mat[flexibles[i]-1]+1)
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

    return template_energy

def get_Z_matrix(pose, optsolution, out_file, resfile, template_energy, flexibles):
    copy_pose = Pose()
    copy_pose.assign(pose)
    score_fxn = create_score_function('talaris2014')
    task_design = TaskFactory.create_packer_task(pose)
    f = open(out_file,'w')
    f.write("Etemp = "+str(template_energy)+"\n")
    task_design.initialize_from_command_line()
    parse_resfile(pose, task_design, resfile)
    pose.update_residue_neighbors()
 #   png = create_packer_graph(pose, score_fxn, task_design)  #Uncomment for latest Pyrosetta versions
    rotsets = RotamerSets()
    ig = pack_rotamers_setup(pose, score_fxn, task_design, rotsets)
#    ig = InteractionGraphFactory.create_and_initialize_two_body_interaction_graph(task_design, rotsets, pose, score_fxn, png)  #Uncomment for latest Pyrosetta versions
    setup_IG_res_res_weights(pose, task_design, rotsets, ig) #Comment for latest Pyrosetta versions
    mat = optsolution
    
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


def compare_seq(seq1,seq2):
  seq_mut=[]
  if seq1==seq2:
    return seq_mut
  else:
    for i in range(0,len(seq1)-1):
      if seq1[i] != seq2[i]:
        seq_mut.append((i,seq2[i]))
  return seq_mut

def read_from_fasta(sequence_file):
  sequences={}
  
  if os.path.exists( os.getcwd() + '/' + sequence_file ) and sequence_file:
    seq_file = open(sequence_file,'r')
    lines = seq_file.readlines()
    read_seq=False
    native=(lines[0][1:-1],lines[1][:-1])
    for line in lines:
      if line[0] == '>':
        read_seq=False
        sequences_name=line[1:-1]
      elif len(line) != 0:
        mut_fasta=line[:-1]
        mutations=compare_seq(native[1],mut_fasta)
        read_seq=True
      if read_seq:
        sequences[sequences_name]=mutations
    return sequences

def mutation_docking(pdb_file, sequence_file, jobs):
  if os.path.exists( os.getcwd() + '/' + pdb_file ) and pdb_file:
    init()
    pose=Pose()
    pose=pose_from_pdb(pdb_file)
    input_file_name=os.getcwd() + '/' + pdb_file.split('.pdb')[0]
    ## create a score function and a scorefile
    scorefxn = create_score_function('talaris2014')
    #score_file_name = input_file_name + "_mut.sc"
    #score_file = open(score_file_name,'w')
    
    ## parse the fasta sequence file into a dictionary : sequences[#position]=[list of mutation]
    sequences={}
    sequences=read_from_fasta(sequence_file)
    
    ## parse the flexible file to get the flexibles residues
    flexibles=open("flexibles.receptor",'r')
    flexibles=flexibles.readlines()[0]
    flexibles=[int(i) for i in flexibles.split()]
    
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
      for mut_tuple in sequences[mut]:
        index = mut_tuple[0] + 1
        aa = mut_tuple[1]
        mutator=MutateResidue(int(index), one_to_three[aa])
        mutator.apply(mut_pose)
        
      mut_folder=os.getcwd()+"/"+mut
      if not os.path.exists(mut_folder):
        os.mkdir(mut_folder)
        
      mut_pose.dump_pdb(mut_folder+"/"+mut+'.pdb')

      #~ pose_packer = standard_packer_task(mut_pose)
      #~ pose_packer.restrict_to_repacking()
      #~ packmover = PackRotamersMover(scorefxn, pose_packer)
      #~ packmover.apply(mut_pose)
        
      ## Minimization on the mutable pose. Useless ?
      minmover.apply(mut_pose)
      
      mut_pose.dump_pdb(mut_folder+"/"+mut+'_min.pdb')

      io = PDBIO()
      pdb = PDBParser().get_structure(mut, mut_folder+"/"+mut+"_min.pdb")
      chain_name=[]
      for chain in pdb.get_chains():
        io.set_structure(chain)
        io.save(mut_folder+'/'+pdb.get_id() + "_" + chain.get_id() + ".pdb")
        chain_name.append(chain.get_id())
      
      ## Separate the two chains ex: E_I or A_B. Need to rename if more than two chains ?
      pose_prot_1=Pose()
      pose_prot_2=Pose()
      pose_prot_1=pose_from_pdb(mut_folder+'/'+pdb.get_id() + "_" + chain_name[0] + ".pdb")
      pose_prot_2=pose_from_pdb(mut_folder+'/'+pdb.get_id() + "_" + chain_name[1] + ".pdb")
      #repack a faire ???? 
      
      ## Minimise the lonely partners
      minmover.apply(pose_prot_1)
      minmover.apply(pose_prot_2)
      
      ## Write the total score of the mutant complex, partner 1 and partner2
      #command=str(mut+' '+str(scorefxn(mut_pose))+' '+str(scorefxn(pose_prot_1))+' '+str(scorefxn(pose_prot_2))+'\n')
      #score_file.write(command)
      
      #### Compute FULL SCP Matrix and 
      compute_interactions(pose_prot_1,'full.resfile', mut_folder+'/receptor.mat')
      command=["Convertor",mut_folder+"/receptor.mat",mut_folder+"/receptor.LG"]
      call(command)
      
      compute_interactions(pose_prot_2,'full.resfile', mut_folder+'/ligand.mat')
      command=["Convertor",mut_folder+"/ligand.mat",mut_folder+"/ligand.LG"]
      call(command)
      
      ###### Calculate the optimal solution and optimal energy.
      command=["toulbar2",mut_folder+"/receptor.LG"]
      tb2out=check_output(command)
      tb2out=tb2out.split('\n')

      for line in tb2out:
        line_split=line.split()
        if ("Optimum:" in line_split) and ("Energy:" in line_split):
            OptEnergy=float(line_split[3])
        elif ("Optimum:" in line_split) :
			OptSolution=line_split[1].split('-')
      OptSolution = [int(i) for i in OptSolution]
      
      #~ mutables=[]
      #~ for tuples in sequences.values():
        #~ for mut in tuples:
            #~ if len(mut)>0:
                #~ mutables.append(mut[0])
      #~ mutables = sorted(list(set(mutables)))
      #~ print mutables
      
      Etemp=get_template_energy(pose_prot_1,OptSolution,OptEnergy,"full.resfile",flexibles)		
      print Etemp
      exit(1)
      #partners=chain_name[0]+'_'+chain_name[1]
      #sample_docking(mut_pose, partners, 1, 1, jobs, mut_folder+"/"+mut)
  else:
    "ERROR: PDB FILE NOT IN THE FOLDER"

        

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
        
mutation_docking(pdb_file,sequence_file, jobs)

