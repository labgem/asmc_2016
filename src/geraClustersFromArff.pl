($headersFileName, $outArffFileName, $msaFileName) = @ARGV;


open (IF, $msaFileName);
$id = "";
while (<IF>){
        chomp $_;
	if (substr($_, 0, 1) eq ">"){
		$id = substr($_,1,10);
	}else{
		$sequence{$id} = $_;
	}
}
close (IF);


open (IF, $headersFileName);
$i = 0;
while (<IF>){
	chomp $_;
	$ids[$i] = $_;
	$i++;	
}
close (IF);


open (IF, $outArffFileName);
$read = 0;
$cluster = 0;
open (OF, ">cwclusters2");
while (<IF>){
	chomp $_;
	if ($read == 1){
		@val = split (/,/, $_);
		$last = $#val;
		$val[$last] =~ s/cluster//;
		print OF $val[$last]." ".$ids[$val[0]]."\n";
#		print ">".$ids[$val[0]]." ".$val[$last]."\n".$sequence{$ids[$val[0]]}."\n";
		print $val[$last]." ".$ids[$val[0]]."\n";
	}else{
		if ($_ eq "\@data"){
			$read = 1;
		}
	}
}
close (IF);
close (OF);
