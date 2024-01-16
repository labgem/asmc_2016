#
# $Id: run_modeller.sh,v 1.1 2011/02/10 $
#
# CMDS: ./run_modeller.sh REPERTORY PERC_identity
# EX:   ./run_modeller.sh ../PF00108 30
# DO: 	Run modeller automatically
# NEED in REPERTORY:
# 	- templates.txt: listing PDBs using as "templates" for MODELLER + chain ID
#			1M4T	A
#			2D3T	A
#
#	- FILEs.pdb: 1M4T.pdb, 2D3T.pdb  
#
# 	- targets.txt: containing all fasta sequences as "targets" for MODELLER
# 	- FILEs.fasta: A6X326.fasta A8ZUE2.fasta



REP=$1
perc_iden=$2
templates="templates.txt"
targets="targets.txt"
SRC=	# <path_to_repository>/asmc_2016/src

N=100  # Number of models in each chunk

echo "........ Treating Directory:" $REP

#
cd $REP || exit 1;

if [ ! -f $templates ]
    then
    echo "$templates not found" >&2
    exit 1
fi


if [ -f  $targets ]
    then
    echo "STARTING ........ Treating repertory:" $REP >&2
    split -l $N $targets $targets.split.
    for i in  $targets.split.*
      do
      if [ ! -d  TMP.$i ]
	  then
	  mkdir TMP.$i
      fi
      (

      cd TMP.$i
      ln -f -s ../$i
      while read j
	do
	ln -f -s ../$j.fasta
	done<$i

      ln -f -s ../$templates
      while read j k
	do
	ln -f -s ../$j.pdb
	done<$templates

      $SRC/modpy.sh python $SRC/model.py $i $templates $perc_iden
      )
    done
else
    echo "$targets not found" >&2
    exit 1
fi

exit 0
