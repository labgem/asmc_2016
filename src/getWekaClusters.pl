#!/usr/bin/perl -w

($arffFileName, $arffOutFileName,  $header) = @ARGV;


open (IF, $header);
$i = 1;
while (<IF>){
    if(/>(\S+)/) {
	$ids[$i] = $1;
	$i++;	
    }
}
close (IF);

$nSeq = $i - 1;
$i=0;
open (IF, $arffFileName);
while (<IF>){
    if ($i >= 1){
	$seq = $_;
	chomp $seq;
	$seq =~ s/,//g;
	$seq[$i]=$seq;
	$i++;
    }	
    if (substr($_, 0, 5) =~ m/data/){
	$i = 1;
    }
}



open (IF, $arffOutFileName);
$index=1;


while (defined ($l = <IF>)){
    if ($l =~ /^(\d+)\s+(\d+)/) {
	$c = $2;
	$index= $1 + 1;
	print ($c,"\t", $ids[$index], "\t", $seq[$index],"\n");
    }
}
close IF;
