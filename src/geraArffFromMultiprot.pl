#!/usr/bin/perl -w 


use strict;
use Getopt::Long;
use vars qw(@ARGV);


my ($msaLength, $msaFileName) = @ARGV;

print "\@relation family\n";

for (my $i=1; $i<=$msaLength; $i++){
	print "\@attribute aa$i {A,R,N,D,C,Q,E,G,H,I,L,K,M,F,P,S,T,W,Y,V,-, X}\n";
}

print "\@data\n";

open (IF, $msaFileName)||die $msaFileName, $!;

my $print = 0;
my @res;
my $i;
while (<IF>){
    chop $_;
    if ($print == 1){
	@res = split (//, $_);
	for ($i=0; $i<$#res; $i++){
	    if ($res[$i] eq "-"){ $res[$i] = '-'; }
	    print $res[$i].",";
	}
	if ($res[$i] eq "-"){ $res[$i] = '-'; }
	print $res[$#res]."\n";
	$print = 0;
    }
    if (substr($_, 0, 1) eq ">"){
	$print = 1;
	}
}
close (IF);
