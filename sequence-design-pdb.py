#!/usr/bin/env python

import os
import optparse
import re
import subprocess
from subprocess import *
import numpy as np

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
    


def rot_trans(pose, partners, flexibles, translation, rotation , jobs, out, mut):

    starting_p = Pose()
    starting_p.assign(pose)
    dock_jump = 1

    setup_foldtree(starting_p, partners, Vector1([1]))

    scorefxn_talaris = create_score_function('talaris2014')
    scorefxn_talaris(starting_p)
    
    
    interface = Interface(1)
    interface.distance(5.0)
    interface.calculate(pose)
    contact = interface.contact_list()
    centroid_muts=[] ## array for the mutables centroids
    centroid_flexs=[] ## array for the flexibles centroids i.e centroid of flexible 1, flexible 2 etc...
    for res in resmuts: ## Calculate the centroid of flexibles of (of res in resnames)
      centroid_flexs.append(centroid([pose.residue(i).xyz("CA") for i in sorted(list(set(contact[int(res)])))]))
      centroid_muts.append(pose.residue(int(res)).xyz("CA"))
    centroid_flex=centroid(centroid_flexs) ## calculate the centroid of flexibles centroids
    centroid_muts=centroid(centroid_muts) ## calculate the centroid of mutables
    X_axis = [centroid_flex[0]-centroid_muts[0],centroid_flex[1]-centroid_muts[1],centroid_flex[2]-centroid_muts[2]] ## Calculate the axis 

    
    ### Define the translation axis by taking the axis between
    ### the centroid of flexibles and the centroid of mutable ?
    axis=rosetta.numeric.xyzVector_Real()
    axis.assign(X_axis[0],X_axis[1],X_axis[2])
    
    trans_pert = RigidBodyTransMover(starting_p,dock_jump)
    trans_pert.trans_axis(axis)
    
    #[2.5134999999999934, -1.6895000000000024, -5.932000000000002]
    
    rot_plus=RigidBodyDeterministicSpinMover()
    rot_plus.spin_axis(axis)
    
    movemap = MoveMap()
    movemap.set_jump(1, True)
    movemap.set_bb(True)
    
    tolerance = 0.01
    min_type = "dfpmin"
    minmover = MinMover(movemap, scorefxn_talaris, min_type, tolerance, True) 
    
    
    minmover.apply(starting_p)
    starting_p.dump_pdb(out+"/PDB/"+mut+"_min.pdb")
    
    compute_interactions(starting_p,'full.resfile', out+"/LG/"+mut+'_min.LG')
    
    command=["toulbar2",out+"/LG/"+mut+'_min.LG']
    tb2out=check_output(command)
    tb2out=tb2out.split('\n')
    for line in tb2out:
        line_split=line.split()
        if ("Optimum:" in line_split) and ("Energy:" in line_split):
            OptEnergy=float(line_split[3])
        elif ("Optimum:" in line_split) :
            OptSolution=line_split[1].split('-')
            OptSolution = [int(i) for i in OptSolution]
    get_Z_matrix(starting_p,OptSolution,OptEnergy,"full.resfile",flexibles,out+"/ZLG/"+mut+'_min.LG')
                                                                           
    
    jd = PyJobDistributor(out+'/PDB/'+mut, jobs, scorefxn_talaris)
    jd.native_pose = starting_p
    counter=1
    ## REGARDER LA BOUCLE A FAIRE
    ## ROTATION ==> [-Teta,Teta] toutes les teta step
    ## Translation ==> [-Delta, Detla] toutes les delta step
    ## Nombre de job en fonction ? ou faire une boucle directement sur le nombre de rotation translation (plus facile)?
    
    Teta=np.range(-rotation*2,rotation*2,0.4)
    Delta=np.range(-translation*2,translation*2,0.4)
    pose = Pose()
    for delta in Delta:
      trans_pert.step_size(delta) # Set the translation size
      for teta in Teta:
        pose.assign(starting_p) # Reload the initial pose
        rot_plus.angle_magnitude(rotation) # Set the translation angle
        trans_pert.apply(pose)
        rot_pert.apply(pose)
        pose.dump_pdb(out+'/PDB/'+mut+'_'+str(counter)+pdb)
        counter += 1
      
    #~ while not jd.job_complete:
      #~ starting_p.assign(jd.native_pose)
      #~ perturb.apply(starting_p)
      #~ 
      #~ compute_interactions(starting_p,'full.resfile', out+"/LG/"+mut+"_"+str(counter)+".LG")
      #~ 
      #~ command=["toulbar2",out+"/LG/"+mut+"_"+str(counter)+".LG"]
      #~ tb2out=check_output(command)
      #~ tb2out=tb2out.split('\n')
      #~ for line in tb2out:
        #~ line_split=line.split()
        #~ if ("Optimum:" in line_split) and ("Energy:" in line_split):
            #~ OptEnergy=float(line_split[3])
        #~ elif ("Optimum:" in line_split) :
            #~ OptSolution=line_split[1].split('-')
      #~ OptSolution = [int(i) for i in OptSolution]
      #~ 
      #~ get_Z_matrix(starting_p, OptSolution,OptEnergy, "full.resfile", flexibles,out+"/ZLG/"+mut+"_"+str(counter)+".LG")
                #~ 
      #~ 
      #~ jd.output_decoy(starting_p)
      #~ counter += 1

def compute_interactions(pose, resfile, out):
    score_fxn = create_score_function('talaris2014')
    task_design = TaskFactory.create_packer_task(pose)
    task_design.initialize_from_command_line()
    parse_resfile(pose, task_design, resfile)
    pose.update_residue_neighbors()
    png = create_packer_graph(pose, score_fxn, task_design)  #Uncomment for latest Pyrosetta versions
    rotsets = RotamerSets()
    ig = pack_rotamers_setup(pose, score_fxn, task_design, rotsets)
    ig = InteractionGraphFactory.create_and_initialize_two_body_interaction_graph(task_design, rotsets, pose, score_fxn, png)  #Uncomment for latest Pyrosetta versions
    g = open(out,'w')
    ener=0.0
    g.write("MARKOV\n")
    g.write(str(ig.get_num_nodes())+'\n')
    domain=StringIO.StringIO()
    scope=StringIO.StringIO()
    unary_terms=StringIO.StringIO()
    binary_terms=StringIO.StringIO()
    number_of_functions=0
    for res1 in range(1,ig.get_num_nodes()+1):
        number_of_functions += 1
        domain_res=str(rotsets.rotamer_set_for_moltenresidue(res1).num_rotamers())
        domain.write(domain_res+' ')
        scope.write("1 "+str(res1-1)+'\n')
        unary_terms.write(domain_res+'\n')
        for i in range(1, rotsets.rotamer_set_for_moltenresidue(res1).num_rotamers()+1):
            resname = rotsets.rotamer_set_for_moltenresidue(res1).rotamer(i).name3()
            ener=ig.get_one_body_energy_for_node_state(res1,i) + reference_energy[resname]
            unary_terms.write(str(-ener)+' ')
        unary_terms.write('\n')
    domain.write('\n')

    for res1 in range(1,ig.get_num_nodes()+1):
        for res2 in range(res1+1,ig.get_num_nodes()+1):
            if (ig.get_edge_exists(res1, res2)):
                number_of_functions += 1
                scope.write("2 "+str(res1-1)+" "+str(res2-1)+"\n")
                domain_res=str(rotsets.rotamer_set_for_moltenresidue(res1).num_rotamers()*rotsets.rotamer_set_for_moltenresidue(res2).num_rotamers())
                binary_terms.write(domain_res+'\n')
                for i in range(1, rotsets.rotamer_set_for_moltenresidue(res1).num_rotamers()+1):
                    nres1=rotsets.rotamer_set_for_moltenresidue(res1).rotamer(i).name3()
                    for j in range(1, rotsets.rotamer_set_for_moltenresidue(res2).num_rotamers()+1):
                        ener=ig.get_two_body_energy_for_edge(res1,res2,i,j)
                        nres2=rotsets.rotamer_set_for_moltenresidue(res2).rotamer(j).name3()
                        binary_terms.write(str(-ener)+' ')
                    binary_terms.write('\n')

    g.write(domain.getvalue())
    g.write(str(number_of_functions)+'\n')
    g.write(scope.getvalue())
    g.write(unary_terms.getvalue())
    g.write(binary_terms.getvalue())
    domain.close()
    scope.close()
    unary_terms.close()
    binary_terms.close()
    g.close()
    
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
        nres=rotsets.rotamer_set_for_moltenresidue(flexibles[i]).rotamer(int(mat[flexibles[i]-1]+1)).name3()
        template_energy-=ig.get_one_body_energy_for_node_state(flexibles[i],int(mat[flexibles[i]-1]+1))+reference_energy[nres]
        for j in range(0,len(mat)):
            if (ig.get_edge_exists(flexibles[i], j+1)):
                if (j+1) in flexibles:
                    if (flexibles[i]<j+1):
                        template_energy-=ig.get_two_body_energy_for_edge(flexibles[i],j+1,int(mat[flexibles[i]-1]+1),int(mat[j]+1))
                else:
                    if (flexibles[i]<j+1):
                        template_energy-=ig.get_two_body_energy_for_edge(flexibles[i],j+1,int(mat[flexibles[i]-1]+1),int(mat[j]+1))
                    else:
                        template_energy-=ig.get_two_body_energy_for_edge(j+1,flexibles[i],int(mat[j]+1),int(mat[flexibles[i]-1]+1))
    
    g = open(out_matrix,'w')
    domain=StringIO.StringIO()
    num_fct=1
    scope=StringIO.StringIO()
    zeroary_terms=StringIO.StringIO()
    zeroary_terms.write("0\n")
    zeroary_terms.write(str(-template_energy)+'\n')
    scope.write('0\n')
    unary_terms=StringIO.StringIO()
    binary_terms=StringIO.StringIO()
    
    for res1 in range(0, len(flexibles)):
        domain_num=rotsets.rotamer_set_for_moltenresidue(flexibles[res1]).num_rotamers()
        domain.write(str(domain_num)+' ')
        scope.write("1 "+str(res1)+'\n')
        unary_terms.write(str(domain_num)+'\n')
        num_fct += 1
        for i in range(1, rotsets.rotamer_set_for_moltenresidue(flexibles[res1]).num_rotamers()+1):
            nres1=rotsets.rotamer_set_for_moltenresidue(flexibles[res1]).rotamer(i).name3()
            unary_ener=ig.get_one_body_energy_for_node_state(flexibles[res1],i) + reference_energy[nres1]
            
            for res2 in range(1,ig.get_num_nodes()+1):
                if res2 not in flexibles:
                    if (ig.get_edge_exists(flexibles[res1], res2)):
                        if (flexibles[res1]<res2):
                            unary_ener+=ig.get_two_body_energy_for_edge(flexibles[res1],res2,i,mat[res2-1]+1)
            unary_terms.write(str(-unary_ener)+' ')
        unary_terms.write('\n')

    binares=[]
    for res1 in range(0, len(flexibles)-1):
        for res2 in range(res1+1, len(flexibles)):
            if (ig.get_edge_exists(flexibles[res1], res2)):
                num_fct += 1
                N= rotsets.rotamer_set_for_moltenresidue(flexibles[res1]).num_rotamers()*rotsets.rotamer_set_for_moltenresidue(res2).num_rotamers()
                scope.write("2 "+str(res1)+' '+str(res2)+'\n')
                binary_terms.write(str(N)+'\n')
                for i in range(1, rotsets.rotamer_set_for_moltenresidue(flexibles[res1]).num_rotamers()+1):
                    for j in range(1, rotsets.rotamer_set_for_moltenresidue(res2).num_rotamers()+1):
                        ener=str(ig.get_two_body_energy_for_edge(flexibles[res1],flexibles[res2],i,j))
                        binary_terms.write(str(ener)+' ')
                    binary_terms.write('\n')



    domain.write('\n')

    binary_terms.close()

    g.write('MARKOV\n')
    g.write(str(len(flexibles))+'\n')
    g.write(domain.getvalue())
    g.write(str(num_fct)+'\n')
    g.write(scope.getvalue())
    g.write(zeroary_terms.getvalue())
    g.write(unary_terms.getvalue())
    domain.close()
    scope.close()
    zeroary_terms.close()
    unary_terms.close()
    g.close()


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
    print sequences
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
      mut_pose=Pose()
      mut_pose.assign(pose)
      print "Processing Mutation:",mut
      resnames=[]
      for mut_tuple in sequences[mut]:
        index = int(mut_tuple[0])
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
        
      if not os.path.exists(mut_folder+'/ZLG'):
        os.mkdir(mut_folder+'/ZLG')
        
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
      compute_interactions(pose_prot_1,'full.resfile', mut_folder+'/LG/receptor.LG')
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
      
      get_Z_matrix(pose_prot_1, OptSolution, OptEnergy, "full.resfile", flexibles_rec, mut_folder+"/ZLG/receptor.LG")	
      
      ##### FOR THE LIGAND
      flexibles_lig_renum=[i-int(chain_length[0]) for i in flexibles_lig]
      compute_interactions(pose_prot_2,'full.resfile', mut_folder+'/LG/ligand.LG')
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
      
      get_Z_matrix(pose_prot_2,OptSolution,OptEnergy,"full.resfile",flexibles_lig_renum,mut_folder+"/ZLG/ligand.LG")		
     
      # FOR THE COMPLEX
      compute_interactions(mut_pose,'full.resfile', mut_folder+"/LG/"+mut+'_min.LG')
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
      
      get_Z_matrix(mut_pose,OptSolution,OptEnergy,"full.resfile",flexibles,mut_folder+"/ZLG/"+mut+'_min.LG')
      
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

