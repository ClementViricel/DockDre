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
    
def centroid(points):
    x_coords = [p[0] for p in points]
    y_coords = [p[1] for p in points]
    z_coords = [p[2] for p in points]
    _len = float(len(points))
    centroid_x = sum(x_coords)/_len
    centroid_y = sum(y_coords)/_len
    centroid_z = sum(z_coords)/_len
    return [centroid_x, centroid_y,centroid_z]

  
def trans_and_rot_to_origin(pdb,center,axis):
  parser = PDBParser()
  structure = parser.get_structure('PDB', pdb)
  chain_list = Selection.unfold_entities(structure, 'C')
  axis=Vector(axis[0],axis[1],axis[2])
  
  r = rotmat(axis,Vector(0,0,1))
  
  for atom in structure.get_atoms():
    coord=atom.get_coord()
    atom.set_coord(map(sub,coord,center))
    
  for atom in structure.get_atoms():
    coord=atom.get_vector()
    new_coord=coord.left_multiply(r)
    new_coord=[i for i in new_coord]
    atom.set_coord(new_coord)
    
  return structure

def Interface_axis(pose, dist, resmuts, score): ## return the interface rotation/translation axis and the center
    
    ## Compute interface at "dist" of distance and store: 
    ## The contact_list=[i := contact of ith residue ] 
    ## The interface_list= [[interface_1],[interface_2]]
    copy_pose=Pose()
    copy_pose.assign(pose)

    score(copy_pose)
    interface = Interface(1)
    interface.distance(dist)
    interface.calculate(copy_pose)
    contact = interface.contact_list()
    interface_residue = interface.pair_list()
    
    centroid_interface=[]
    for interface in interface_residue: ## Compute the axis by taking the whole interfaces.
      centroid_interface.append(centroid([copy_pose.residue(res).xyz('CA') for res in interface])) ## store the centroids of the residue in the interfaces.
    interface_axis= map(sub,centroid_interface[0],centroid_interface[1])
    center=centroid_interface[0]
    return (interface_axis, center)

def rot_trans(pose, partners, flexibles, translation, rotation , trans_step, rot_step, out, mut, resmuts,scorefxn):

    copy_pose = Pose()
    copy_pose.assign(pose)
    dock_jump = 1

    scorefxn(copy_pose)

    compute_interactions(copy_pose,'full.resfile', out+"/LG/"+mut+'_min.LG',scorefxn)
    copy_pose.dump_pdb(out+"/PDB/"+mut+"_min.pdb")
    os.rename(out+"/LG/"+mut+'_min.uai', out+"/UAI/"+mut+'_min.uai')
    
    interface_axis_center=Interface_axis(copy_pose, 10, resmuts, scorefxn)
    interface_axis=interface_axis_center[0]
    center=interface_axis_center[1]
    
    command=["toulbar2",out+"/LG/"+mut+'_min.LG',"-w="+out+"/SOL/"+mut+'_min.sol']
    tb2out=call(command)
    
    io = PDBIO()
    structure=trans_and_rot_to_origin(out+"/PDB/"+mut+"_min.pdb", center, interface_axis) ## Translate the complex to have the center of rotation on origin
    Teta=np.arange(-rotation,rotation+rot_step,rot_step)
    Delta=np.arange(-translation,translation+trans_step,trans_step)
    counter=1
    out_command=open("matrix-command.txt",'ab')
    three_dim=[Vector(0,0,1),Vector(0,1,0),Vector(1,0,0)] ## Three dimentional vector, forgot ones ?
    for delta in Delta: ## loop for trans
      for teta in Teta: ## loop for rotations along the three axis
		  for i in three_dim:
			structure_copy=structure.copy()
			rotation=rotaxis(np.pi/180*teta, i) ## Matrice de rotation de vecteur (0,0,1) et de centre (0,0) (vecteur "interface_axis" et centre "center")
			translation=np.array((0, 0, delta), 'f') ## translation
			chain_list = Selection.unfold_entities(structure_copy, 'C') #Selectionne les chaines A et B (il faut bien les renommer avant)
			chain_list[1].transform(rotation, translation) # Rotate et Translate la chaine B seulement
			io.set_structure(structure_copy)
			io.save(out+'/PDB/'+mut+'_'+str(counter)+"_("+str(delta)+")T_("+str(teta)+")R.pdb")
			command="python command_d_t.py --out "+out+" --count "+str(counter)+" --delta "+str(delta)+" --teta "+str(teta)+" --mut "+mut+"\n"
			out_command.write(command)
			counter += 1
    out_command.close()
        

def compute_interactions(pose, resfile, out, score_fxn):
    copy_pose=Pose()
    copy_pose.assign(pose)
   
    ## Minimization
    movemap = MoveMap()
    movemap.set_jump(1, True)
    movemap.set_bb(True) ## Move Backbone
    tolerance = 0.001 ## Set the tolerance, minimum will be in tolerance of a local min
    min_type = "dfpmin" ## Quasi Newton
    minmover = MinMover(movemap, score_fxn, min_type, tolerance, True) 
    minmover.apply(copy_pose)
    #relax=FastRelax(score_fxn)
    #relax.apply(copy_pose)
    
    task_design = TaskFactory.create_packer_task(copy_pose)
    task_design.initialize_from_command_line()
    parse_resfile(copy_pose, task_design, resfile)
    copy_pose.update_residue_neighbors()
    png = create_packer_graph(copy_pose, score_fxn, task_design)
    rotsets = RotamerSets()
    ig = pack_rotamers_setup(copy_pose, score_fxn, task_design, rotsets)
    ig = InteractionGraphFactory.create_and_initialize_two_body_interaction_graph(task_design, rotsets, copy_pose, score_fxn, png) 
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
    
def mutation_rot_trans(pdb_file, seq_file, translation_size, rotation_size, translation_step, rotation_step,beta):
  if os.path.exists( os.getcwd() + '/' + pdb_file ) and pdb_file:
    if os.path.exists ("matrix-command.txt"):
      os.remove("matrix-command.txt")
    if beta==True:
        init('-ex1 false -ex1aro false -ex2 false -ex2aro false -beta')
        scorefxn = create_score_function('beta_july15')
    elif beta==False:
        init('-ex1 false -ex1aro false -ex2 false -ex2aro false')
        scorefxn = create_score_function('talaris2014')

    pose=Pose()
    pose=pose_from_pdb(pdb_file)
    setup_foldtree(pose,"A_B",Vector1([1]))
    input_file_name=os.getcwd() + '/' + pdb_file.split('.pdb')[0]

    ## parse the fasta sequence file into a dictionary : sequences[#position]=[list of mutation]
    sequences={}
    sequences[pdb_file.split('.pdb')[0]]=[]
    if os.path.exists( os.getcwd() + '/' + seq_file ) and seq_file:
        r = re.compile("([0-9]+)([a-zA-Z]+)")
        seqfile=open(seq_file,'r')
        for seq in seqfile.readlines():
            sequences[seq[:-1]]=[r.match(i).groups() for i in seq.split('_')]
    print sequences
    ## parse the flexible file to get the flexibles residues
    flexibles_rec=open("flexibles.receptor",'r')
    flexibles_rec=flexibles_rec.readlines()[0]
    flexibles_rec=[int(i) for i in flexibles_rec.split()]
    
    flexibles_lig=open("flexibles.ligand",'r')
    flexibles_lig=flexibles_lig.readlines()[0]
    flexibles_lig=[int(i) for i in flexibles_lig.split()]
    
    flexibles=sorted(flexibles_rec+flexibles_lig)
    
    ## First minimisation
    #movemap = MoveMap()
    #movemap.set_jump(1, True)
    #movemap.set_bb(True)
    #tolerance = 0.001
    #min_type = "dfpmin"
    #minmover = MinMover(movemap, scorefxn, min_type, tolerance, True) 
    #minmover.apply(pose)
    #pose.dump_pdb(os.getcwd()+"/"+input_file_name+"_min_0.pdb")
    
    
    ###### Mutation Loop #### Can be paralelle ?
    for mut in sequences.keys():
      mut_pose=Pose()
      mut_pose.assign(pose)
      print "Processing Mutation:",mut
      resmuts=[]
      for mut_tuple in sequences[mut]:
        index = int(mut_tuple[0])
        aa = mut_tuple[1]
        resmuts.append(index)
        mutator=MutateResidue(int(index), one_to_three[aa])
        mutator.apply(mut_pose)
        
      mut_folder=os.getcwd()+"/"+mut
      if not os.path.exists(mut_folder):
        os.mkdir(mut_folder)
        
      if not os.path.exists(mut_folder+'/PDB'):
        os.mkdir(mut_folder+'/PDB')
        
      if not os.path.exists(mut_folder+'/SOL'):
        os.mkdir(mut_folder+'/SOL')

      if not os.path.exists(mut_folder+'/UAI'):
        os.mkdir(mut_folder+'/UAI')

      if not os.path.exists(mut_folder+'/LG'):
        os.mkdir(mut_folder+'/LG')
      
      if os.path.exists(mut_folder+'/score.sc'):
        os.remove(mut_folder+'/score.sc')
        
      mut_pose.dump_pdb(mut_folder+"/PDB/"+mut+'.pdb')
      
      io = PDBIO()
      pdb = PDBParser().get_structure(mut, mut_folder+"/PDB/"+mut+".pdb")
      chain_name=[]
      chain_length=[]
      for chain in pdb.get_chains():
        io.set_structure(chain)
        io.save(mut_folder+'/PDB/'+pdb.get_id() + "_" + chain.get_id() + ".pdb")
        chain_length.append(len(chain))
        chain_name.append(chain.get_id())
      partners=chain_name[0]+'_'+chain_name[1]
      
      ## Separate the two chains. Need to rename in A and B
      pose_prot_1=Pose()
      pose_prot_2=Pose()
      pose_prot_1=pose_from_pdb(mut_folder+'/PDB/'+pdb.get_id() + "_" + chain_name[0] + ".pdb")
      pose_prot_2=pose_from_pdb(mut_folder+'/PDB/'+pdb.get_id() + "_" + chain_name[1] + ".pdb")
            
      ###### Compute FULL SCP matrix, Calculate the optimal solution and optimal energy.
      ##### FOR THE RECEPTOR
      compute_interactions(pose_prot_1,'full.resfile', mut_folder+'/LG/receptor.LG',scorefxn)
      os.rename(mut_folder+'/LG/receptor.uai', mut_folder+'/UAI/receptor.uai')
      command=["toulbar2",mut_folder+"/LG/receptor.LG","-w="+mut_folder+"/SOL/receptor.sol"]
      tb2out=call(command)
      
      ##### FOR THE LIGAND
      flexibles_lig_renum=[i-int(chain_length[0]) for i in flexibles_lig]
      compute_interactions(pose_prot_2,'full.resfile', mut_folder+'/LG/ligand.LG',scorefxn)
      os.rename(mut_folder+'/LG/ligand.uai', mut_folder+'/UAI/ligand.uai')
      command=["toulbar2",mut_folder+"/LG/ligand.LG","-w="+mut_folder+"/SOL/ligand.sol"]
      tb2out=call(command)

      ## Loop for trans rot
      rot_trans(mut_pose, partners, flexibles, translation_size, rotation_size, translation_step, rotation_step, mut_folder,mut,resmuts,scorefxn)
    
      print "Finish Processing Mutation:",mut
  else:
    "ERROR: PDB FILE NOT EXISTING"

        

parser=optparse.OptionParser()
parser.add_option('--pdb', dest = 'pdb_file',
    default = '',    
    help = 'Protein comple in PDB format' )

parser.add_option('--seq', dest = 'seq_file',
    default = '',    
    help = 'Sequences to map' )

parser.add_option( '--trans', dest='translation_size' ,
    default = 2.0,    
    help = 'Size of translation segment')
  
parser.add_option( '--rot', dest='rotation_size' ,
    default = 3.0,   
    help = 'Size of rotation segment')
    
parser.add_option( '--trans_step', dest='translation_step' ,
    default = 0.4,   
    help = 'Size of translation steps')
    
parser.add_option( '--rot_step', dest='rotation_step' ,
    default = 2.0,   
    help = 'Size of rotation steps')

parser.add_option( '--beta', dest='beta' ,
    default = False,   
    help = 'Beta scoring function')       
(options,args) = parser.parse_args()

pdb_file=options.pdb_file
sequence_file = options.seq_file
translation_size=float(options.translation_size)
rotation_size=float(options.rotation_size)
translation_step=float(options.translation_step)
rotation_step=float(options.rotation_step)
beta=options.beta


################# MUTATION, PDB and MATRIX PRODUCTION #############

start_time=datetime.now()
print "Score function beta:",beta
mutation_rot_trans(pdb_file, sequence_file, translation_size, rotation_size, translation_step, rotation_step,beta)
end_time=datetime.now()-start_time
print end_time
