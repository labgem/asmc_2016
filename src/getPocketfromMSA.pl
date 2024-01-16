#!/usr/bin/perl -w 



use strict;
use Getopt::Long;
use vars qw(@ARGV);




my $usage = "\nUSAGE:
 getPocketfromMSA -ref <REFPDBID> <MSA> <positionsFile>

    ";

my $aide = 0;
my $ref=undef;
unless (GetOptions("help", \$aide,
		   "ref=s", \$ref)){
    die $usage;
}
if(@ARGV != 2) {
    warn $usage;
    exit(1);
}


my $MSAFile = shift @ARGV;
my $positionsFile =  shift @ARGV;


if ($aide) {
    print "$usage";
    exit(1);
}
if(!(defined($ref))){
    warn "Argument -ref is mandatory\n$usage";
    exit(1);
}


sub loadMSA($) {
    my $file = shift;
    my %tab = ();
    if(!(open(I,$file))){
	warn "$file : $!";
	return undef;
    }
    while(defined(my $l=<I>)){
	if($l=~/^(\S+)\s+(\S+)/) {
	    if(defined( $tab{$1})) {
		$tab{$1} = $tab{$1} . $2;
	    } else {
		$tab{$1} = $2;
	    }
	}
    }
    close I;
    return(\%tab);
}

sub  loadPositions($) {
    my $file = shift;
    my @tab = ();
    if(!(open(I,$file))){
	warn "$file : $!";
	return undef;
    }
    while(defined(my $l=<I>)){
	if($l=~/^(\d+)$/){
	    push(@tab,$1);
	} else {
	    warn "$file : not a valid line : $l";
	}
    }
    close I;
    return(\@tab);
}


sub computeCoord($) {
    my $seq = shift;
    my @tab;
    my $x;
    my $l;
    my $c;
    my $index=0;
    for($x=0;$x<length($seq)-1;$x++){
	$c = substr($seq,$x,1);
	if (!($c eq "-")) {
	    $tab[$index] = $x;
	    $index++;
	}
    }
    return(\@tab);
}

my $MSA = loadMSA($MSAFile);
my $positions = loadPositions($positionsFile);
my $coordonnees;

if(defined($MSA->{$ref})) {
    $coordonnees=computeCoord($MSA->{$ref});
} else {
    warn "$ref not found in $MSAFile";
    exit(1);
}


for my $id (keys(%$MSA)) {
    print ">$id\n";
    for my $pos (sort {$a <=> $b} @$positions) {
	print substr($MSA->{$id},$coordonnees->[$pos-1],1);
    }
    print "\n";
}
