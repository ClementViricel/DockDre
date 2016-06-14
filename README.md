# DockDre

Script that mutate the residue of --pdb.
Mutations are taken in the sequence file --seq (the file must be written in raw)
A folder will be create for each mutation with the mutation name.
Hierarchy of the folder is:

                              MUT
                               |
                      PDB_____SOL ____ LG ____(ZLG useless for now)

Example Command line:

python sequence-design.py --pdb 1PPF.pdb --seq 1PPF.test.seq

python Jay-Z --pdb 1PPF.pdb --seq 1PPF.test.seq

It will take 5 hours to run this two lines. (The first command is quick and the second is long)

ENJOY !
