use strict;

use Math::Complex;

my $usage = "$0 <subsitution matrix_file> <AlignClusterfile> <AlignOtherClusters>";

my $matrixFileName =  shift @ARGV;
my $clusterFileName = shift @ARGV;
my $othersFileName =  shift @ARGV;
my $sizeSequence = 0;
my @globalAlignment = ();
my $numClusters = 0;
my $numOthers = 0;
my @LLi;
my @LLRi;

# Read WAG substitution
my %mm;

my %vue;

readMatrixFile($matrixFileName,\%mm);


#Input data
open (IF, $clusterFileName)|| die "$! : $usage";
while (<IF>){
    my $j = 1;    
    if (/^>(\S+)/) {
	$vue{$1}=1; #exclude entry from others clusters.
    } else {
	s/\s//g;
	chomp $_;
	my @seq = split (//, $_);	
	$numClusters++;
	$sizeSequence = length $_;	
	foreach my $aa (@seq){
	    $globalAlignment[$numClusters][$j] = $aa;
	    $j++;
	}
    }
}

close (IF);
open (IF, $othersFileName);
my $todo;

my $x = $numClusters;
while (<IF>){
    my $j = 1;
    if (/^>(\S+)/) {
	if($vue{$1}) {
	    $todo = 0;
	} else {
	    $todo = 1;
	}
    } elsif ($todo) {
	s/\s//g;
	chomp $_;
	my @seq = split (//, $_);
	$numOthers++;
	$x++;
	foreach my $aa (@seq){
	    $globalAlignment[$x][$j] = $aa;
	    $j++;
	}
    }
}

close (IF);

sub computeLLi ($$) {
    my $alignment = shift;
    my $j =shift;
    my $A ;
    my $B ; 
    my $C ; 
    my $D ; 
    my $AA ; 

    my %count;
    my %countGlobal;
    my $aa_count_cluster = 0 ; # number of non gap position in the column
    my $aa_count_other = 0;    # number of non gap position in the column
    my $LLi;


    
    for(my $i=1;$i<=$numClusters;$i++) {
	if (not(defined($alignment->[$i][$j]))) {
	    warn $i," ",$alignment->[$i][$j];
	} else {
	    $count{$alignment->[$i][$j]}++;
	}
    }	
    for(my $i=1+$numClusters;$i<=$numClusters+$numOthers;$i++) {
	$countGlobal{$alignment->[$i][$j]}++;
    }
    
#Find the most frequent residues
    my @aas =  ('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y');

    my $max = 0; 
    my $maxaa = "";

    foreach my $aa (@aas){
	$aa_count_cluster += $count{$aa};
	$aa_count_other   +=  $countGlobal{$aa};
	if ($count{$aa} > $max){
	    $max = $count{$aa};
	    $AA = $aa;
	}
    }
  
  
    $A  = $max;
    $B = $countGlobal{$AA};
    #$C = $numClusters; # count only regular aa (no gap)
    $C =  $aa_count_cluster;
    #$D = $numOthers;
    $D =  $aa_count_other ;



    if(($A==0)||(($C+$D)==0)) {
	return 0;
    }
    
    my $e0 = ($A + $B) / ($C + $D);
    my $e1 = $C * $e0 ;
    my $e2 = $D * $e0 ;
    
    
    if ($B != 0){
	$LLi = 2 * (($A * ln($A/$e1)) + ($B * ln($B/$e2)));
    }else{
	$LLi = 2 * ($A * ln($A/$e1));
    }
    return($LLi);
}



for(my $j = 1; $j<=$sizeSequence;$j++) {
    $LLi[$j] = computeLLi(\@globalAlignment,$j);
}



#Permutations
my $number = $numClusters + $numOthers;
for (my $j=1; $j<=$sizeSequence; $j++){ 
    my $sum = 0;
    for (my $k=1; $k<=100; $k++){
	my @bkp = ();
	for (my $i=1; $i<=$number; $i++){
	    $bkp[$i][$j] = $globalAlignment[rand($number)+1][$j];
	}
	$sum += computeLLi(\@bkp,$j);
    }
    $sum /= 100;
    
    if($sum == 0) {
	warn "sum nulle";
	$LLRi[$j] = $LLi[$j] / 0.0001;
    } else {
	$LLRi[$j] = $LLi[$j] / ($sum);
    }
}

for (my $j=1; $j<=$sizeSequence; $j++){
    my $sum = 0;
    my $dp;
    my @LLisim;

    for (my $k=1; $k<=1000; $k++){
	my @bkp = ();
	for (my $i=1; $i<=$number; $i++){
	    my $r = int (rand ($number)+1);
	    my $r2 = int(rand(10000)) + 1;	    
	    $bkp[$i][$j]= $mm{$globalAlignment[$r][$j].$r2};
	    if(not(defined($bkp[$i][$j]))) {
		$bkp[$i][$j]= "-";
	    }
	}
	$LLisim[$j][$k] = computeLLi(\@bkp,$j);
	$sum+=$LLisim[$j][$k];
    }
    
    $sum /= 1000; 
    for (my $k=1; $k<=1000; $k++){
	$dp += ($LLisim[$j][$k]-$sum)**2;
    }
    
    if($dp == 0) {
	warn "variance nulle";
	$dp =+ 0.0001;
    }
    $dp = sqrt ($dp);

    my $zScore = abs($LLRi[$j] - $sum) / ($dp); 
    print "$j ".$LLRi[$j]." ".$sum." ".$dp." ".$zScore."\n";
}



sub readMatrixFile($$) {
    my $file = shift;
    my $matrix = shift;
    my @aa;
    my @line;
    my @val;
    my $lineNumber= 0;

    open(I,$file);
    while (<I>) {
	chomp;
	next if (/^>/);
	if (/\#/) {
	    @line = split (/\s+/,$_);
	    for (my $i=0;$i<@line; $i++)  {
		push(@aa,$line[$i]);
		$val[$i] = 0;
	    }
	} else {
	    @line = split (/\s+/,$_);
	    $lineNumber++;
	    for (my $i=1;$i<@line; $i++)  {
		my $val = int(100 * $line[$i]);
		for (my $x=$val[$i];$x<=$val[$i]+$val;$x++){
		    $mm{$aa[$i] . $x} = $aa[$lineNumber];
		}
		$val[$i]+= $val;
	    }
	}
    }
}
