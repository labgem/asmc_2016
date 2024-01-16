from sys import *
from Bio.Align import AlignInfo
from Bio import AlignIO


# USAGE: python conserv.py all.cluster 40

#
#
# Write maximal conservation for each position in conservation.dat
#
#

REP=str(argv[3])
allcluster=str(argv[1])

alignment = AlignIO.read(open(allcluster), "fasta")

alig_len = alignment.get_alignment_length()


if len(argv) == 4:
    threshold=int(argv[2])
else:
    threshold=40.0

ref_seq = alignment.get_seq_by_num(0)

cons = open(REP+"/consensus.txt","w")
summary_align = AlignInfo.SummaryInfo(alignment)
consensus = summary_align.dumb_consensus()

cons.write(str(consensus))

my_pssm = summary_align.pos_specific_score_matrix(consensus,chars_to_ignore = ['X'])
max = len(alignment)

# ----------------------
cpsfile= open(REP+"/cps.dat","w")
output = open(REP+"/conservation.dat","w")
print "All conservation % are now written in conservation.dat file"
print "For information: conservation above %d" % threshold, "%:"
for pos in xrange(alig_len):
    perc_max = 0
    letter_max = "-"
    for letter in my_pssm[pos].keys():
        percent = (my_pssm[pos][letter] / max) * 100.0
        if percent > perc_max and letter != "-" :
           perc_max = percent
           letter_max = letter
        if percent > threshold and letter != "-":
            print "%d %s %3.2f%s" % (pos+1, letter, percent, ' %')
            cps = "%s\n" % (pos+1)
            cpsfile.write(cps)
    conservation = "%d %s %3.2f%s\n" % (pos+1, letter_max, perc_max, ' %') 
    output.write(conservation)
