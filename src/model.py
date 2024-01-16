# 
#
#	DATE: May 4th 2010
#   USAGE: modpy.sh python model.py
#
#
#
#

import sys
import os
import glob

from modeller import *
from modeller.automodel import *
from modeller.scripts import complete_pdb
from modeller.parallel import *
#from Bio import SeqIO

# Use 8 CPUs in a parallel job on this machine
#j = job()
#j.append(local_slave())
#j.append(local_slave())
#j.append(local_slave())
#j.append(local_slave())
#j.append(local_slave())
#j.append(local_slave())
#j.append(local_slave())
#j.append(local_slave())

# ARGUMENTS:

# List of fasta files
targets=str(sys.argv[1])
templates=str(sys.argv[2])
perc_identity=float(sys.argv[3])


# Cleaning ....
for filename in glob.glob('xrays_ali.*'):
   os.remove (filename)

# Writing the outputs:
output = open('models.dat',"a")


 
log.verbose() # request verbose output
env = environ() # create a new MODELLER environment to build this model in
env.io.atom_files_directory = './:../atom_files/' # read topology
env.libs.topology.read(file='$(LIB)/top_heav.lib')  # read parameters

# Need :
# - templates.txt: list of Xray PDBS files
# - targets.txt: list of fasta sequences 

# Cleaning the repertory in case of....
for filename in glob.glob('xrays_ali.*'):
   os.remove (filename)

# Writing the outputs:
output = open('models.dat',"a")

CODE = []
CHAIN = []
for line in file(templates):
    line = line.split()
    CODE.append(line[0])
    CHAIN.append(line[1])



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
               dendrogram_file='xrays.tree',
               alignment_type='tree', # If 'progresive', the tree is not
                                      # computed and all structues will be
                                      # aligned sequentially to the first
               feature_weights=weights, # For a multiple sequence alignment only
                                        # the first feature needs to be non-zero
               improve_alignment=True, fit=True, write_fit=write_fit,
               write_whole_pdb=whole, output='ALIGNMENT QUALITY')

aln.write(file='xrays_ali.pap', alignment_format='PAP')
aln.write(file='xrays_ali.ali', alignment_format='PIR')

aln.salign(rms_cutoff=1.0, normalize_pp_scores=False,
           rr_file='$(LIB)/as1.sim.mat', overhang=30,
           gap_penalties_1d=(-450, -50), gap_penalties_3d=(0, 3),
           gap_gap_score=0, gap_residue_score=0, dendrogram_file='1is3A.tree',
           alignment_type='progressive', feature_weights=[0]*6,
           improve_alignment=False, fit=False, write_fit=True,
           write_whole_pdb=False, output='QUALITY')

KNOWNS = []
for seq in aln:
   KNOWNS.append(seq.code)

# PART II: align target sequence on multiple structures

TARGETS = []
for line in file(targets):
    line = line.split()
    TARGETS.append(line[0])

for target in TARGETS:
   print target
   a = alignment(env, file=target+'.fasta', alignment_format='FASTA')
   a.write(file='target.ali', alignment_format='PIR')
   # Read aligned structure(s):
   aln = alignment(env)
   aln.append(file='xrays_ali.ali', align_codes='all')
   aln_block = len(aln)

   # Read aligned sequence(s):
   aln.append(file='target.ali', align_codes=target)

   # Structure sensitive variable gap penalty sequence-sequence alignment:
   aln.salign(output='', max_gap_length=20,
           gap_function=True,   # to use structure-dependent gap penalty
           alignment_type='PAIRWISE', align_block=aln_block,
           feature_weights=(1., 0., 0., 0., 0., 0.), overhang=0,
           gap_penalties_1d=(-450, 0),
           gap_penalties_2d=(0.35, 1.2, 0.9, 1.2, 0.6, 8.6, 1.2, 0., 0.),
           similarity_flag=True)
   aln.write(file='target-mult.ali', alignment_format='PIR')
   aln.write(file='target-mult.pap', alignment_format='PAP')
   
   percent=0.00
   nb = len(aln)
   for seq in aln: # sequence identity
       print seq.get_sequence_identity(aln[nb-1])
       if seq.get_sequence_identity(aln[nb-1]) > percent and seq.get_sequence_identity(aln[nb-1])!= 100.00:
           percent=seq.get_sequence_identity(aln[nb-1])
   print 'percent', percent
      
   if percent > perc_identity:
       a = automodel(env, 
                     alnfile='target-mult.ali',
                     knowns=(KNOWNS), 
                     sequence=target,
                     assess_methods=(assess.DOPE, assess.GA341)
                     )
       a.starting_model = 1
       a.ending_model = 5
       a.make()
       ok_models = filter(lambda x: x['failure'] is None, a.outputs)

       # DOPE function
       ok_models.sort(lambda a,b: cmp(a['DOPE score'], b['DOPE score']))
       mini = ok_models[0]

       # CLEANING

       os.rename(mini['name'], target + ".pdb")
       nb = len(aln)
       output.write(target)

       for filename in glob.glob(target + '.V*'):
          os.remove(filename)
       for filename in glob.glob(target + '.D*'):
          os.remove(filename)
       for filename in glob.glob(target + '.B*'):
          os.remove(filename)
       for filename in glob.glob(target + '.{r,s}*'):
          os.remove(filename)
       os.remove (target+".ini")
       os.remove (target+".rsr")
       os.remove (target+".sch")

   for seq in aln: # sequence identity
       ident_seq = "   %s :  %3.2f  " % (seq.code, seq.get_sequence_identity(aln[nb-1]))
       output.write(ident_seq)
   output.write("\n")

os.remove("target-mult.ali")
os.remove("target-mult.pap")
os.remove("target.ali")
if os.path.isfile("xrays.tree"):
   os.remove("xrays.tree")
os.remove("xrays_ali.ali")
os.remove("xrays_ali.pap")


   
