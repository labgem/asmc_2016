from sys import *
from Bio.Align import AlignInfo
from Bio import AlignIO


# USAGE: python conserv_pos.py all.cluster 40

#
#
# Write maximal conservation for each position in conservation.dat
#
#

# REP=str(argv[3])
allcluster=str(argv[1])

alignment = AlignIO.read(open(allcluster), "fasta")

alig_len = alignment.get_alignment_length()


if len(argv) == 3:
    threshold=int(argv[2])
else:
    threshold=40.0

summary_align = AlignInfo.SummaryInfo(alignment)
consensus = summary_align.dumb_consensus()

my_pssm = summary_align.pos_specific_score_matrix(consensus,chars_to_ignore = ['X'])
max = len(alignment)

# ----------------------

pre_pos = 0
pos_dublicate = 'false'
for pos in xrange(alig_len) :
	perc_max = 0
	letter_max = "-"
	for letter in my_pssm[pos].keys():
		percent = (my_pssm[pos][letter] / max) * 100.0
		if percent > perc_max and letter != "-" :
			perc_max = percent
			letter_max = letter
			if percent > threshold and letter != "-" :
				# print "%d %s %3.2f%s" % (pos+1, letter, percent, ' %')##############
				if pre_pos == pos+1 :
					pos_dublicate = 'true'
					break
				pre_pos = pos+1
	if pos_dublicate == 'true' :
		break

print pos_dublicate
