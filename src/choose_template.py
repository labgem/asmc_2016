import Bio.PDB
import numpy
import os
from Bio.PDB.PDBParser import PDBParser
import sys
import math
# USAGE: python choose_template.py REP yes/no
# yes: take the structure in complex / no: take the best resolution structure no in complex


REP=str(sys.argv[1])
option = str(sys.argv[2])
list_templates = open(REP+"/"+"templates.txt", "rU")

templates=[]
chain=[]
for line in list_templates:
   line = line.split()
   templates.append(line[0])
   chain.append(line[1])

parser=PDBParser(PERMISSIVE=0)

min_resolution=10 
min_template="TOTO"
templates_complex=[]

for template in templates:
   structure = parser.get_structure(template, REP+"/"+template+".pdb")
   keywords=structure.header['keywords']
   journal=structure.header['journal']
   if(keywords.rfind("complex")==-1 and journal.rfind("COMPLEX")==-1 and option=="no" ):
      resolution=structure.header['resolution']
      if(resolution<min_resolution):
         min_resolution=resolution
         min_template=template
   else:
      templates_complex.append(template) # if keyword == 1, means structure in complex

# Select struc with best resolution and beeing in complex (templates_complex)

if(min_template=="TOTO" or option=="yes"):
   for template in templates_complex:
      structure = parser.get_structure(template, REP+"/"+template+".pdb")
      resolution = structure.header['resolution']
      if(resolution<min_resolution):
         min_resolution=resolution
         min_template=template


index=templates.index(min_template)
min_chain=chain[index]

print min_template,"  ",min_chain
