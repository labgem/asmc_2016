#!/usr/bin/perl -w 


use strict;
use Getopt::Long;
use vars qw(@ARGV);




my $usage = "\nUSAGE:
 maxEntropie.pl <fastafile>";


my ($fastaFile) = @ARGV;

open (IF, $fastaFile);

my @mat;
my $l;

while(defined($l=<IF>)){
    next if ($l =~ /^>(\S+)/);
    chomp $l;
    my @t = split('',$l);
    push (@mat,\@t);
}
close IF;

my $netat = 21;
my $maxEntropie = 0;
my $i=0;
my $n=0;
my @comptes;

for my $tab (@mat) {
    $n++;
    $i=0;
    if(!(defined($comptes[$i]))){
	$comptes[$i] = {};
    }
    for my $aa (@{$tab}){
	if(!(defined( $comptes[$i]->{$aa} ))) {
	    $comptes[$i]->{$aa} = 1;
	} else {
	    $comptes[$i]->{$aa}++;
	}
	$i++;
    }
}

for (my $col=0;$col<$i;$col++) {
    my $H=0;
    for my $aa (keys %{$comptes[$col]}) {
	my $aa_value = $comptes[$col]->{$aa};
	$H -=  $aa_value/$n * log($aa_value / $n) / log(2);
    }
    $H = log($netat)/log(2) - $H;
    if ($H > $maxEntropie){
	$maxEntropie = $H;
    }
}

printf("%.2f\n",$maxEntropie);

0;
