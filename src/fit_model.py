# 
#
#	DATE: Karine May 4th 2010
#   USAGE: modpy.sh python model.py
#
#
#
#


import sys
import os
import glob

sys.path.append('') # <path_to_repository>/asmc_2016/modeller/<path_to_python>
# e.g: sys.path.append('/home/eelisee/asmc/v0/modeller/linux-x86_64-generic/lib/x86_64-intel8/python2.7')

from modeller import *
from modeller.automodel import *
from modeller.scripts import complete_pdb
#from Bio import SeqIO

log.verbose() # request verbose output
env = environ() # create a new MODELLER environment to build this model in
env.io.atom_files_directory = './:../atom_files/' # read topology
env.libs.topology.read(file='$(LIB)/top_heav.lib')  # read parameters

# Need :
# - template.txt: list of Xray PDBS files

# ARGUMENTS:

# Repertory to work on it:
REP = str(sys.argv[1])
PDBREF = str(sys.argv[2])
CHAINREF= str(sys.argv[3])


print "...Working Directory.................", REP


# Cleaning the reperory in case of....
##for filename in glob.glob(REP + 'xrays_ali.*'):
##   os.remove (filename)


CODE = []
CHAIN = []

for line in file(REP+"/templates.txt"):
    line = line.split()
    CODE.append(REP+line[0])
    CHAIN.append(line[1])

CODE.append(REP+"/"+PDBREF)
CHAIN.append(CHAINREF)

# PART I: salign
# Illustrates the SALIGN multiple structure/sequence alignment
# Align structures

aln = alignment(env)

for i in range(len(CODE)): 
    code=CODE[i]
    chain=CHAIN[i]
    mdl = model(env, file=code, model_segment=('FIRST:'+chain, 'LAST:'+chain))
    aln.append_model(mdl, atom_files=code, align_codes=code+chain)

for (weights, write_fit, whole) in (((1., 0., 0., 0., 1., 0.), False, True),
                                    ((1., 0.5, 1., 1., 1., 0.), False, True),
                                    ((1., 1., 1., 1., 1., 0.), True, False)):
    aln.salign(rms_cutoff=3.5, normalize_pp_scores=False,
               rr_file='$(LIB)/as1.sim.mat', overhang=30,
               gap_penalties_1d=(-450, -50),
               gap_penalties_3d=(0, 3), gap_gap_score=0, gap_residue_score=0,
               dendrogram_file=REP+'xrays.tree',
               alignment_type='tree', # If 'progresive', the tree is not
                                      # computed and all structues will be
                                      # aligned sequentially to the first
               feature_weights=weights, # For a multiple sequence alignment only
                                        # the first feature needs to be non-zero
               improve_alignment=True, fit=True, write_fit=write_fit,
               write_whole_pdb=whole, output='ALIGNMENT QUALITY')

aln.write(file=REP+'all_pdbs_ali.pap', alignment_format='PAP')
aln.write(file=REP+'all_pdbs_ali.ali', alignment_format='PIR')

aln.salign(rms_cutoff=1.0, normalize_pp_scores=False,
           rr_file='$(LIB)/as1.sim.mat', overhang=30,
           gap_penalties_1d=(-450, -50), gap_penalties_3d=(0, 3),
           gap_gap_score=0, gap_residue_score=0, dendrogram_file='1is3A.tree',
           alignment_type='progressive', feature_weights=[0]*6,
           improve_alignment=False, fit=False, write_fit=True,
           write_whole_pdb=False, output='QUALITY')

