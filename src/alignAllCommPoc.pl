#!/usr/bin/perl -w 

# 	$Id: alignAllCommPoc.pl,v 1.2 2011/02/09 09:35:39 artigue Exp $	

use strict;
use Getopt::Long;
use vars qw(@ARGV);

my $usage = "\nUSAGE:
 alignAllCommPoc.pl  -of <Fasta outFile1> -op <Position outFile2> <templateFileName> <resPocFileName> <listFileName> <delFile>

    ";

my $aide = 0;
my $ref=undef;
my $saida=undef;
my $matches=undef;


my %del;
my %poc;

unless (GetOptions("help", \$aide,
		   "of=s",\$saida,
		   "op=s",\$matches)){
    die $usage;
}

if ($aide) {
    print "$usage";
    exit(1);
}

if(!(defined($saida))||(!defined($matches))){
    warn "both of and op are mandatories\n$usage";
    exit(1);
}

my ($templateFileName, $resPocFileName, $listFileName, $delFile) = @ARGV;

sub trim($){
    my $string = shift;
    $string =~ s/^\s+//;
    $string =~ s/\s+$//;
    return $string;
}

if(defined($delFile)) {
    open (IF, $delFile);
    while (<IF>){
        chomp $_;
        $del{$_} = 1;
    }
    close (IF);
}

open (IF, $resPocFileName)||die "$resPocFileName : $!";
while (<IF>){
	chomp $_;
	$poc{$_} = 1;
}
close (IF);

open (IF, $listFileName) ||die "$listFileName : $!";
open (OF, ">$saida");
open (OF2, ">$matches");

while (my $fileName = <IF>){
    chomp $fileName;
    if(!(open (IF2, "$fileName")))  {
	warn "$fileName : $!";
	next;
    }
    my %poc2 = (); 
    my %poc3 = ();
    my $ok = 0;

    while (my $line = <IF2>){
	chomp $line;
	if ((substr($line, 0, 5) eq "Match")&&($ok == 0)){
	    $ok = 1;
	}elsif (substr($line, 0, 3) eq "End"){
	    $ok = 2;
	}elsif ($ok == 1){
	    my $resa = trim(substr($line, 0, 8));
	    my $resb = trim(substr($line, 9, 8));
	    my ($chain1, $res1, $resNum1) = split (/\./, $resa);
	    my ($chain2, $res2, $resNum2) = split (/\./, $resb);
	    if (defined($poc{$resNum1})){
		$poc2{$resNum1} = $res2;
		$poc3{$resNum1} = $resNum2;
	    }
	}
    }
    close (IF2);
    my $seq = "";
    my $i = 1;
    my $ACNUM=$fileName;
    $ACNUM =~ s/.*\///;
    $ACNUM =~ s/\..*//;
    
    foreach my $key (sort {$a <=> $b} keys %poc){
	if ((exists $poc2{$key})&&(! exists $del{$i})){
			print OF2 "$ACNUM $i ".$poc3{$key}."\n";
			$seq .= $poc2{$key};
		    }else{
			$seq .= "-";
		    }
	$i++;
    }
    print OF ">$ACNUM\n$seq\n";
}
close (IF); 
close (OF);
close (OF2);


0;
