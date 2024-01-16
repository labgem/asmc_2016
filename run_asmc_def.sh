#!/usr/bin/env bash


ROOT=	# <path_to_repository>/asmc_2016
BIN=$ROOT/bin
SRC=$ROOT/src

export PYTHONPATH=$BIN/CoreBio-0.5.0/:$PYTHONPATH
weblogo=$BIN/weblogo-3.0/weblogo
wekajar=$BIN/weka-3-4-13/weka.jar
multiprot=$ROOT/MultiProtInstall/multiprot.Linux
multiprotParams=$ROOT/MultiProtInstall/params.txt
fpocket=$BIN/fpocket

matrix=$BIN/HAson.mat
cutoff=0.01 # WEKA CUTOFF

##############################################
#  ASMC	Active Site Clustering and Modeling  #
#					     #
#  	April 2016			     #
##############################################
	
pg=`basename $0`

usage() {
    echo 
    echo
    echo "USAGE: $pg [ -d DIR ] -models FILENAME -pdbref FILENAME -chain CHAIN_ID -pockets N"
    echo
    echo "-d DIR             :  Working Directory to store results (default ./)"
    echo
    echo "-models FILENAME   :  List of homology models files (Format PDB)"
    echo 
    echo "-pdbref            : PDB file used as a reference structure. The best reference structure is"
    echo "                     the one with a ligand inside the active site or a cavity. Otherwise, pick up "
    echo "                     the structure with the best resolution or one with no loops,"
    echo "                     side-chains or fragments missing."
    echo 
    echo
    echo "-chain CHAIN_ID    : Chain for the reference structure."
    echo
    echo "-pockets N         :  The number of pockets you want to analyse from 1 to N. fpocket"
    echo "                      ranks detected cavities from the biggest to the smallest one." 
    echo "	                Generally, the active pocket corresponds to the biggest one or "
    echo "                      the more conserved one. The user should check that the"
    echo "                      pocket entered corresponds to the active site."
    echo "                     (see fpocket manual for more information)"
    echo
    echo "OUTPUTS:"
    echo
    echo "DIR/Multiprot (Multiprot output)"
    echo "   Results of structural pairwise alignment of each structures with the reference structure"
    echo
    echo "DIR/1REF_out (1REF of the reference structure"
    echo "   Results of fpocket program)"
    echo
    echo "DIR/Pocket_[0-N]"
    echo "   list_aa_pocket.dat          : List of amino acids from PDBREF that belong to POCKET"
    echo "   align.pos                   : Correspondance between numerotation of the pocket and"
    echo "                                 numerotation of amino acids in the sequence"
    echo "   align.fasta                 : Fasta format of all sequences"
    echo "   consensus.txt               : Consensus sequence of POCKET"
    echo "   conservation.dat            : Conservation for each position of the pocket."
    echo "   clusters_after_fusion.fasta : Fasta format of all sequences belonging to one cluster"
    echo "   tree.nw                     : Tree of the most clusters (newick) could be read by: iTol (http://itol.embl.de/)"
    echo "   tree.nexus                  : Tree of the most clusters (nexus)"
    echo "   tree_all_cluster.nw         : Tree of all the clusters (newick)"
    echo "   allseq.png                  : Logo sequence of the pocket for all sequences"
    echo "   xx.cluster                  : Fasta format of sequences belonging to the xx.cluster"
    echo "   xx.png                      : Logo sequence of the pocket for the xx cluster"
    echo "   xx_pvalue.dat               : pvalues for each position fo the pocket"
    echo
    echo
    echo
    echo "   Example: $pg -d /tmp -models models.txt -pdbref 2VDJ.pdb -chain A -pockets 3" 
    echo
    echo 
    echo
    echo 
    echo "   Extraction and structural alignment of pockets from homology models"
    echo "   Generate ASMC tree"
    echo "   Analysis of SDPs & CPs (Melo-Minardi et al., Bioinformatics. 2010 Dec 15;26(24):3075-82)"
    echo 
}


log() {
    echo "$1" >&2
}

die(){
	usage >&2
	exit 1
}


REP=
fModels=
PDBREF=
CHAIN=
POCKET=

while :
  do
  case x"$1" in
      x-d)    shift; REP="$1"
	  shift
	  ;;
      x-models) shift; fModels="$1"
	  shift
	  ;;
      x-pdbref) shift; PDBREF="$1"
	  shift
	  ;;
      x-chain) shift; CHAIN="$1"
	  shift
	  ;;
      x-pockets)shift;  POCKET="$1"
	  shift
	       ;;
	x) break
		;;
	x*) die
		;;
	esac
done

case x"$REP" in
    x) REP="./"
	;;
esac

case x"$fModels" in
    x) echo "-models Mandatory" >&2
	die
	;;
esac
if [ ! -f $fModels ] 
    then 
    echo "Not a regular file: $fModels" >&2
    die
fi

if [ `cat $fModels | wc -l` -lt 500 ] ; then
	cluster_size=5
elif [ `cat $fModels | wc -l` -lt 1000 ]; then
	cluster_size=10
else 
	cluster_size=15
fi


case x"$PDBREF" in
    x) echo "-pdbref Mandatory" >&2
	die
	;;
esac

if [ $PDBREF = "pick" ]
then
   PDBREF_CHAIN=`python $SRC/choose_template.py $REP yes`
   PDBREF=`echo $PDBREF_CHAIN | awk '{print $1}'`
else
   PDBREF=`echo $PDBREF | sed 's/.pdb$//'`
fi

if [ ! -f "$PDBREF.pdb" ]
    then
    echo "$PDBREF.pdb not found" &>2
    die
fi



case x"$CHAIN" in
    x) echo "-chain Mandatory" >&2
    die
    ;;
esac

if [ "$CHAIN" = "pick" ]
   then
   CHAIN=`echo $PDBREF_CHAIN | awk '{print $2}'`
fi 


if [ ! -e "$REP"/"$PDBREF"\_fit.pdb ] ; then echo ".... Creating _fit " ; $SRC/modpy.sh python $SRC/fit_model.py $REP $PDBREF.pdb $CHAIN ; fi
mv  $REP/$PDBREF\_fit.pdb $REP/\1REF.pdb


case x"$POCKET" in
    x) echo "-pockets Mandatory" >&2
	die
	;;
esac

TMP=/tmp/asmc-$$
mkdir $TMP
trap "rm -rf $TMP" 0

MultiProt=$REP/MultiProt
saida="align.fasta"
matches="align.pos"


log `date +"%m-%d-%Y"`

log
log "ASMC analysis  START"
log
log "REP = $REP"
log "Models = $fModels"
log "PDBREF = $PDBREF CHAIN = $CHAIN"
log "POCKETS = $POCKET"
log


if [ ! -d "$REP" ] 
    then
    log "Creating $REP"
    log
    mkdir $REP
fi

log "Multiprot Analysis : START"
log

if [ ! -d "$MultiProt" ] 
    then
    log "Creating $MultiProt"
    log
    mkdir $MultiProt
fi


ModelsRoot=`dirname $fModels`
# PDBID=`basename \1REF.pdb`


while read i
  do
  R=`dirname $i`
  j=`basename $i`

  x=`echo $i | cut -c 1,1`
  if [ "$x" == "/" ] 
      then
      ORI=$R
  else
      ORI=$ModelsRoot/$R
  fi
 
  if [ -f "$MultiProt"/"$i".sol.res ] 
      then
      log "MultiProt : Skip already computed $i"
  else
      log "MultiProt : Compute $i"
      cp \1REF.pdb   $TMP/
      cp $ORI/$j $TMP/
      (
	  cd  $TMP
	  $multiprot \1REF.pdb $j
      )
      mv $TMP/2_sol.res $MultiProt/$j.sol.res
  fi
  done< $fModels

log "Multiprot Analysis END"
log
log "Fpocket Analysis START"
log


# EXTRACTING POCKET FROM ref.pdb
if [ -d  "$REP"/\1REF_out ]
    then
    log  "fPOCKET : Skip already computed pockets for 1REF"
else
    log "fPOCKET : Compute  1REF"
    if [ ! -f "$REP"/\1REF.pdb ]
	then
	cp \1REF.pdb $REP/
	trap "rm $REP/\1REF.pdb" 0
    fi
    (
	cd $REP
	$fpocket  -f \1REF.pdb
    )
fi


log "fPOCKET Analysis END"
log

i=0
while [ $i -lt $POCKET ]
  do
  if [ -f "$REP"/\1REF_out/pockets/pocket"${i}"\_atm.pdb ] ; then
	  log  "POCKET $i START"
	  DEST=$REP/Pocket_$i
	  [ -d $DEST ] || mkdir $DEST
	  [ -d $DEST/LOG ] || mkdir $DEST/LOG
	  
          cp list_aa_pocket.dat $DEST/ 
	  ls $REP/MultiProt/*.sol.res > $TMP/MPList_$i
	  $SRC/alignAllCommPoc.pl -of $DEST/$saida -op $DEST/$matches $PDBREF.pdb $DEST/list_aa_pocket.dat $TMP/MPList_$i 2>$DEST/LOG/alignAllCommPoc.e
	  (
		  cd $DEST
		  nb_aa=`wc -l list_aa_pocket.dat | awk '{print $1}'`

		  $SRC/geraArffFromMultiprot.pl $nb_aa $saida > msa_pocket.arff 2> LOG/geraArffFromMultiprot.e

		  log "WEKA : START"  

		  for nb_run in `seq 1 1` ;
			do
				if [ -f  msa_pocket_arff.tree ] && [ -f  msa_pocket_arff.out ] && [ -f  Weka_clusters ] ; then
				  log "Skipped"
				  break
				else 
				  log "Run WEKA"
				  java -cp $wekajar -Xmx500M weka.clusterers.Cobweb -t msa_pocket.arff -A 1.0 -C $cutoff      > msa_pocket_arff.tree 
				  java -cp $wekajar -Xmx500M weka.clusterers.Cobweb -t msa_pocket.arff -A 1.0 -C $cutoff -p 1 > msa_pocket_arff.out
				  perl $SRC/getWekaClusters.pl msa_pocket.arff  msa_pocket_arff.out $saida | sort -k1n >  Weka_clusters
				fi
				clusterWeka=`awk '{print $1}' Weka_clusters | sort -n | uniq`

				for c in $clusterWeka
				 do
					awk '($1 == c) {print ">" $2 "\n" $3}' c=$c Weka_clusters > $c.fasta
					CONSERV_POS=`python $SRC/conserv_pos.py $c.fasta 20`
					rm $c.fasta
					if [ "$nb_run" -lt 1 -a "$CONSERV_POS" = "true" ] ; then
						cutoff=`bc -l <<<"$cutoff/2"`
						rm msa_pocket_arff.tree msa_pocket_arff.out Weka_clusters #OUTPUTs of run_weka
						break
					fi
				 done
				if [ "$nb_run" -lt 5 -a "$CONSERV_POS" = "true" ] ; then
					continue
				fi
				break
			done 
		  rm msa_pocket.arff
		  log "WEKA : END"

		  log "Merge Clusters : START"
		  
		  perl $SRC/merge_clusters_on_maxGain.pl msa_pocket_arff.tree Weka_clusters > clusters_after_fusion  2> LOG/clusters_after_fusion.err
	          perl $SRC/merge_clusters_on_maxGain.pl -newick msa_pocket_arff.tree clusters_after_fusion > tree.nw 2> tree_all_cluster.err
		  perl $SRC/merge_clusters_on_maxGain.pl -nexus msa_pocket_arff.tree clusters_after_fusion > tree.nx 2> tree.err

                  # Generate the full ASMC tree - with leaves
                  list_group=`awk '{print $1}' clusters_after_fusion | sort -n | uniq`
                  cp tree.nw tree-0.nw
                  l=0
                  k=0

                  for i in $list_group; do l=$(( $k+1 )); list=`awk 'BEGIN {printf "%s", "("} $1=='$i'{printf "%s%1s",$2,","} END {printf "%s", ")"}' clusters_after_fusion` ; j="$list"; echo $j ; sed 's/('$i')/'$j'/g' tree-$k.nw | sed 's/,)/)/g' | sed 's/)(/),(/g' > tree-$l.nw ; k=$(( $k+1 )); done

                  mv tree-$l.nw tree_full.nw
                  if [ -f tree-1.nw ]; then
                  rm tree-*.nw
                  fi

		 log "Merge Clusters : END"
		 log

		 awk '{print ">" $2 "\n" $3}' clusters_after_fusion >  clusters_after_fusion.fasta 
						  
		 clusterList=`awk '{print $1}' clusters_after_fusion | sort -n | uniq -c | awk '($1 >= n) {print $2}' n=$cluster_size`
		 N=`echo $clusterList | wc -w`
		 log "Number of Clusters $N"
		 log


		  # FORMATING 
		  if [ "$N" -gt 0 ]
		  then
		  echo -n > clusters_after_fusion_and_$cluster_size.fasta
		  echo -n > clusters_after_fusion_and_$cluster_size
		  for c in $clusterList
			do
			awk '($1 == c) {print ">" $2 "\n" $3}' c=$c clusters_after_fusion >>  clusters_after_fusion_and_$cluster_size.fasta
			awk '($1 == c) {print $0}'             c=$c clusters_after_fusion >>  clusters_after_fusion_and_$cluster_size
		  done  
			  
			  # MAKING THE TREE OF MOST POPULATED CLUSTERS
		  perl $SRC/merge_clusters_on_maxGain.pl -newick msa_pocket_arff.tree clusters_after_fusion_and_$cluster_size >  tree_$cluster_size.nw    2> LOG/tree_$cluster_size.e
		  perl $SRC/merge_clusters_on_maxGain.pl -nexus  msa_pocket_arff.tree clusters_after_fusion_and_$cluster_size  > tree_$cluster_size.nexus 2> LOG/tree_$cluster_size.e

			  #CALCUL CONSENSUS SEQUENCE OF THE POCKET AND % OF CONSERVATION 
		  python $SRC/conserv.py  clusters_after_fusion.fasta 75 ./
		  
		  
			  # GENERATE LOGO FOR ALL SEQUENCES 
		  nb_seq=`wc -l  clusters_after_fusion | awk '{print $1}'`
		  webLogoY=`$SRC/maxEntropie.pl  clusters_after_fusion.fasta `
		  $weblogo -f  clusters_after_fusion.fasta -o allseq.png -t "All" -P "$nb_seq sequences " -F png -A protein --resolution 300 --color-scheme chemistry --errorbars no --composition none -S $webLogoY
		  convert allseq.png allseq.ppm
			  
		  echo  >  sdps.dat
		  for c in $clusterList
			do
			if [ -e pvalue.dat ] ; then rm pvalue.dat ; fi
			awk '($1 == c) {print ">" $2 "\n" $3}' c=$c clusters_after_fusion > $c.cluster
			nb_seq=`awk '($1 == c)'  c=$c  clusters_after_fusion | wc -l | awk '{print $1}'`
			webLogoY=`$SRC/maxEntropie.pl  $c.cluster`
			$weblogo -f $c.cluster -o $c.png -t "Cluster $c" -P "$nb_seq sequences " -F png -A protein --resolution 300 --color-scheme chemistry --errorbars no --composition none -S $webLogoY
			convert $c.png $c.ppm
			perl $SRC/loglikelihood.pl $matrix $c.cluster clusters_after_fusion.fasta > $c.zscore 2>LOG/loglikelihood$c.e
			rm $c.cluster
			R --vanilla --slave --args $c < $SRC/graph.R
		  done
		  R --vanilla --slave --args tree_$cluster_size.nexus < $SRC/make-asmc-tree.R
		  rm -r *ppm clusters_after_fusion_and_$cluster_size.fasta 
		  fi
		  rm clusters_after_fusion.fasta 
	  )
	fi
  i=`expr $i + 1`
done

exit 0
