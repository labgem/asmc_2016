#!/usr/bin/perl -w 


use strict;
use Getopt::Long;
use vars qw(@ARGV);




my $usage = "\nUSAGE:
 merge_clusters_on_maxGain.pl [-tree | -nexus] arff.tree    arff.clusters

 Analyses WEKA Cobwek outputs (options tree and clusters composition outputs)
 such that nodes are only considered to be split only when at least gain on one
 MSA position  is higher than gains of the upper levels.

outputs : <Node ID> <Seq ID> <Seq pattern>

[-tree] : output information for all nodes of the tree
    <Level> <Node ID> <Parent ID> <Seq ID>,<Seq pattern>
    
    [-nexus] : output the selected tree in nexus format
    
    [-newick] : output the selected tree in newick format" ;

my $newick = 0;
my $nexus = 0;
my $tree = 0;
my $aide = 0;


unless (GetOptions("help", \$aide,
                   "nexus", \$nexus,
                   "newick", \$newick,
		   "tree", \$tree)){
    die $usage;
}

if ($aide) {
    print "$usage";
    exit(1);
}

if($nexus&&$tree) {
   warn "Options tree and nexus are exclusives\n";
   die $usage;

exit(1);
 }

if( @ARGV != 2) {
   warn "2 arguments recquired\n";
   die $usage;
}

my ($arffFileNameTree, $clustersCompFile) = @ARGV;

open (IF, $arffFileNameTree) || die  $arffFileNameTree," : ", $!;

my $netat = 21;
my $parent = -1;

my %nodes;
my %id;
my %origine;
my %childs;
my %mat;
my %clusters;
my %GAIN;
my %EFFECTIF;
my %SPLIT;

 my $l;

while(defined ($l=<IF>)) {
    my $node;
    my $effectif;

    next if ($l =~ /^\s*$/);
    if ($l =~ /[node|leaf]\s*(\d+)\s*\[(\d+)\]/) {
	$node = $1;
	$effectif = $2;
	if(!(defined($nodes{$node}))) {
	    $nodes{$node} = $effectif;
	    $origine{$node} = $parent;
	    $childs{$parent}{$node}=1;
	} 
	$parent = $node;
    }
}
close IF;



open(IF,$clustersCompFile);
while(defined($l=<IF>)){
        next if ($l =~ /^\s*$/);
        if ($l =~ /^\s*(\d+)\s+(\S+)\s+(\S+)/){
	my @t = split(//,$3);
	push (@{$mat{$1}},\@t);
	push (@{$id{$1}},$2);
	$clusters{$1}=1;
    } else {
	die "KO","-- $l\n";
    }
}
close IF;

sub delete_empty_leaves  {
    my ($node) = @_;

    my @children = (keys %{$childs{$node}});

    if(@children > 0){
	for my $nodeF (keys %{$childs{$node}}) {
	    delete_empty_leaves($nodeF);
	}
    } 
    
    @children = (keys %{$childs{$node}});

    if((@children == 0)&&(!(defined($mat{$node})))){
	warn "Deleting empty leave $node\n";
	delete($childs{$origine{$node}}->{$node}) ;
	delete($mat{$node});
	delete($childs{$node});
	delete($nodes{$node});
	delete($origine{$node});
	delete($clusters{$node});
    }
    if((@children > 0)&&(defined ($mat{$node}))){
	my $index = 0;
	for my $p (@{$mat{$node}}) {
	    #push (@{$mat{$child}},$p);
	    #push (@{$id{$child}},$id{$node}[$index]);
	    warn "Deleting  node misplacement ",$id{$node}[$index]," from cluster ",$node,"\n";
#	    warn "Correcting misplacement ",$id{$node}[$index]," from cluster ",$node," to ", $child,"\n";
	    $index++;
	}
	delete $mat{$node};
	delete $id {$node};
    }
}


sub population_nodes {
    my ($node) = @_; 
    for my $nodeF (keys %{$childs{$node}}) {
	population_nodes($nodeF);
	#if(!(defined(@{$mat{$nodeF}}))){
	if(! @{$mat{$nodeF}}){
          die "KO : delete again $nodeF";
	} else {
	    my $index = 0;
	    for my $p (@{$mat{$nodeF}}) {
		push(@{$mat{$node}},$p);
		push(@{$id{$node}},$id{$nodeF}[$index]);
		$index++;
	    }
	}
    }   
}

delete_empty_leaves(0);
population_nodes(0);



my %ENTROPIE;

for my $node (sort {$a <=> $b} keys %nodes) {
    my @comptes=();
    my $n = 0;
    my @n=();
    for my $tab (@{$mat{$node}}) {
	$n++;
	my $i=0;
	if(!(defined($comptes[$i]))){
	    $comptes[$i] = {};
	}

	for my $aa (@{$tab}){
	    if(!(defined( $comptes[$i]->{$aa} ))) {
		$comptes[$i]->{$aa} = 1;
	    } else {
		$comptes[$i]->{$aa}++;
	    }
	    $n[$i]++;
	    $i++;
	}
    }
    for (my $col=0;$col<@comptes;$col++) {
	my $H=0;

	for my $aa (keys %{$comptes[$col]}) {
	    if ($aa eq "-") {
	    } else {
		my $aa_value = $comptes[$col]->{$aa};
		$H -= $aa_value / $n[$col] * log ($aa_value / $n[$col]) / log($netat);
	    }
	}
	$ENTROPIE{$node}[$col]=$H;
    }
#    if ($n != $nodes{$node}) {
#	warn "KO, $node : expected : ",$nodes{$node}, " count ", $n;
#    }
    
    $EFFECTIF{$node}=$n;
}

sub compute_Gains {
    my ($node)=@_;
    my @children = (keys %{$childs{$node}});
    my $split = 0;
    for (my $col=0;$col< @{$ENTROPIE{$node}};$col++) {
 	my $h = $ENTROPIE{$node}[$col];
 	my $effectif = $EFFECTIF{$node};
 	my $max = $GAIN{$node}[$col];
 	for my $nodeF (keys %{$childs{$node}}) {
 	    my $hnext = $EFFECTIF{$nodeF}/$effectif * $ENTROPIE{$nodeF}[$col];
# 	    my $hnext = $ENTROPIE{$nodeF}[$col];
 	    if($h - $hnext > $max){
		$max = $h - $hnext;
 		$split = 1;
		$GAIN{$nodeF}[$col] = $h - $hnext;
	    } else {	    
		$GAIN{$nodeF}[$col] = $max;
	    }
	}
    }
#	for my $nodeF (keys %{$childs{$node}}) {
#	    $GAIN{$nodeF}[$col] = $max;
#	}
    
    
    
    if($split){ 
	$SPLIT{$node}=1;
	for my $nodeF (keys %{$childs{$node}}) {
	    compute_Gains($nodeF);
	}
    } else {
	$SPLIT{$node}=0;
    }
}


for (my $col=0;$col<@{$ENTROPIE{0}};$col++) {
    $GAIN{0}[$col]=0;
}
compute_Gains(0);


sub print_tree  {
    my ($node) = shift;
    my $level = shift;
    my $index=0;
    for my $p (@{$mat{$node}}) {
	print $level,"\t", $node,"\t";
	print $origine{$node},"\t";
	print $id{$node}[$index],"\t";
	print @$p,"\n";
	$index++;
    }
    if($SPLIT{$node}!=0) {
	$level+=1;
	for my $nodeF (keys %{$childs{$node}}) {
	    print_tree($nodeF, $level);
	}
    }
}

sub print_newick  {
    my ($node) = shift;
    print "(",$node;
    for my $nodeF (keys %{$childs{$node}}) {
	print_newick($nodeF);
    }
    print ")";
}


sub print_nexus {
    my $node =shift;
    my $nbclus = scalar keys %{$childs{$node}};
    if($nbclus>1) {
	print help_nexus($node,"");
    } elsif ($nbclus==1) {
	print "(", help_nexus($node,"") ,")";
    } else {
	print "(", 0,")";
    }
}


sub help_nexus  {
    my $node = shift;
    my $chain = shift;
    my $nbclus = scalar keys %{$childs{$node}};
    my $nodeF;
    my $i;
    if($nbclus>1) {
	$chain = ")" . $node . $chain ;
	my @nodes =  (keys %{$childs{$node}});
	for ($i=0;$i<@nodes-1;$i++){
	    $nodeF=$nodes[$i];
	    $chain = "," . help_nexus($nodeF,$chain);
	}
	$nodeF=$nodes[$i];
	$chain = "(" . help_nexus($nodeF,$chain);
    } elsif ($nbclus==1) {
	my @nodes =  (keys %{$childs{$node}});
	$nodeF=$nodes[0];
	$chain = help_nexus($nodeF,$chain);
    } elsif($nbclus==0) {
	$chain = $node . $chain ;
    }
    return($chain);
}

sub print_nodes  {
    my ($node) = @_; 
    if($SPLIT{$node}==0) {
	my $index=0;
	for my $p (@{$mat{$node}}) {
	    print $node,"\t";
	    print $id{$node}[$index],"\t";
	    print @$p,"\n";
	    $index++;
	}
    } else {
	for my $nodeF (keys %{$childs{$node}}) {
	    print_nodes($nodeF);
	}
    }
}



if($tree) {
    print_tree(0,0);
} elsif ($newick) {
    print_newick(0);
    print "\n";
} elsif ($nexus)  {
print "#NEXUS \n"; 
print " \n";
print "Begin trees; \n";
print "\tTree tree1= ";
print_nexus(0);
print ";\nEnd;\n";
} 
else {
    print_nodes(0);
}

0;
