#!/usr/bin/perl


#Perl Script to disect results of vcftools sliding windows and calculate Fold Change and Foldchange proportion of 2 vcftools ouputs. FE teosinte vs Landraces.
#it also take into aconunt files of cordintates to look into particular regions.
# input vcftools --window-pi --window-pi-step
# Author: Miguel Vallebueno  CC BY-NC-SA 4.0
# Date: 2024-08-23





$in="Candidate_genes.txt";
print "Loading candidate list <$in>....\n";
# Open the file or die with an error message
open($fh, '<', $in) or die "Could not open file '$in' $!";





$cmin=100000000000000000000;
$cmax=1;


while (my $line = <$fh>) {
    chomp $line; 
	 
	 @col= split("\t", $line);
	 
	 for($i=$col[3];$i<=$col[4];$i++){ 
		if($col[0] eq "c"){$cand{$col[2]}{$i}=$col[1];
			if($i<$cmin){$cmin=$i}
			if($i>$cmin){$cmax=$i}
		}
		if($col[0] eq "n"){$neig{$col[2]}{$i}=$col[1]}
		
			#print "xxx <$col[1]>   <$i>"; <STDIN>;
	 }
}

# Close the file handle
close($fh);


#print"aaa  <$cmax>  <$cmin>"; <STDIN>;

#$cmax=121494532;
#$cmin=41758978;

$cmax=$cmax+30000;
$cmin=$cmin-30000;

#print"bbb  <$cmax>  <$cmin>";<STDIN>;


$in="Hufford_Domestication_Candidates_chr5.txt";
print "Loading Hufford Domestication loci<$in>....\n";
# Open the file or die with an error message
open($fh, '<', $in) or die "Could not open file '$in' $!";

while (my $line = <$fh>) {
    chomp $line; 
	 
	 @col= split("\t", $line);
	 if($col[1] ne "candidate"){next;}
	 for($i=$col[3];$i<=$col[4];$i++){ 
			$huff{$col[2]}{$i}=$col[0];
	 }
}

# Close the file handle
close($fh);





### loading Outgroup pop Vcftools Diversity estimator outfile

$in="merged_flt_c5.imputed.Teo.Wpi_100.Spi_50.windowed.pi";
print "Loading outgroup vcftols windowed file <$in>....\n";
open($fh, '<', $in) or die "Could not open file '$in' $!";
$mino=1000;
while (my $line = <$fh>) {
	if($line =~ /^CHR/){next;}
    chomp $line; 
	@col= split("\t", $line);
	$key="$col[0]&$col[1]&$col[2]";
	$outv{$key}=$col[4];
	$ALL{$key}=0;
	if($col[4]<$mino){$mino=$col[4]}	
	
}
# Close the file handle
close($fh);


### loading test pop Vcftools Diversity estimator outfile

$in="merged_flt_c5.imputed.LR.Wpi_100.Spi_50.windowed.pi";
print "Loading testpop vcftols windowed file <$in>....\n";
open($fh, '<', $in) or die "Could not open file '$in' $!";
$mint=1000;
while (my $line = <$fh>) {
	if($line =~ /^CHR/){next;}
    chomp $line; 
	@col= split("\t", $line);
	$key="$col[0]&$col[1]&$col[2]";
	$testv{$key}=$col[4];
	$ALL{$key}=0;
	if($col[4]<$mint){$mint=$col[4]}
    	
	#print "<$col[0]> <$col[1]>  <$col[4]>";<STDIN>;
}
# Close the file handle
close($fh);

$mint=$mint-0.000000000001;
$mino=$mino-0.000000000001;

print "Calculating FC....\n";

	#foreach my $key (sort { $a cmp $b } keys %outv) {
	foreach my $key (sort { $a <=> $b }  keys %ALL) {
	
	@col= split(/&/, $key);
	$chr=$col[0];
	$str=$col[1];
	$end=$col[2];	
	
	if(!exists $testv{$key}){$testv{$key}=$mint}
	if(!exists $outv{$key}){$outv{$key}=$mino}
	
	#$fc=$testv{$key}/$outv{$key};
	
	#print "<$key>   <$str> <$cmin> <$cmax>"; <STDIN>;
	
	if($str > $cmin and $end < $cmax ){$fpREG{$key}=$testv{$key}/($outv{$key}+$testv{$key});}
	
	#$fc{$key}=$outv{$key}/$testv{$key};
	if(exists $cand{$chr}{$str} | exists $cand{$chr}{$end}){$fcCA{$key}=$testv{$key}/$outv{$key}; $fpCA{$key}=$testv{$key}/($outv{$key}+$testv{$key}); next;}
	if(exists $neig{$chr}{$str} | exists $neig{$chr}{$end}){$fcNE{$key}=$testv{$key}/$outv{$key}; $fpNE{$key}=$testv{$key}/($outv{$key}+$testv{$key});  next;}
	if(exists $huff{$chr}{$str} | exists $huff{$chr}{$end}){$fcHUF{$key}=$testv{$key}/$outv{$key}; $fpHUF{$key}=$testv{$key}/($outv{$key}+$testv{$key});  next;}
	else{$fcNU{$key}=$testv{$key}/$outv{$key}; $fpNU{$key}=$testv{$key}/($outv{$key}+$testv{$key});}
	
}




$out="Hufford_dom_genes_FC.txt";
open(OUT, ">$out") or die "Could not open file '$out': $!";
print "Printing outputs HUFF....\n";
print OUT "Region\tfc\tfp\n";
foreach my $key (keys %fcHUF){
	@col= split("&", $key);
	$chr=$col[0];
	$str=$col[1];
	$end=$col[2];
	print OUT " $huff{$chr}{$str}\t$fcHUF{$key}\t$fpHUF{$key}\n";
}
close (OUT);



$out="Candidate_FC.txt";
open(OUT, ">$out") or die "Could not open file '$out': $!";
print "Printing outputs CAND....\n";
print OUT "Region\tfc\tfp\n";
foreach my $key (keys %fcCA){
	@col= split("&", $key);
	$chr=$col[0];
	$str=$col[1];
	$end=$col[2];
	print OUT " $cand{$chr}{$str}\t$fcCA{$key}\t$fpCA{$key}\n";
}
close (OUT);

$out="Neigbors_FC.txt";
open(OUT, ">$out") or die "Could not open file '$out': $!";
print OUT "Region\tfc\tfp\n";
print "Printing outputs NE....\n";
foreach my $key (keys %fcNE){
	@col= split("&", $key);
	$chr=$col[0];
	$str=$col[1];
	$end=$col[2];
	print OUT " $neig{$chr}{$str}\t$fcNE{$key}\t$fpNE{$key}\n";
}
close (OUT);

$out="Null_FC.txt";
open(OUT, ">$out") or die "Could not open file '$out': $!";
print OUT "Region\tfc\tfp\n";
print "Printing outputs NULL....\n";
foreach my $key (keys %fcNU){
	print OUT "$key\t$fcNU{$key}\t$fpNU{$key}\n";
}
close (OUT);


$out="REG_FP.txt";
open(OUT, ">$out") or die "Could not open file '$out': $!";
print "Printing outputs REG....\n";
foreach my $key (sort  { $a <=> $b }    keys %fpREG){
	@col= split("&", $key);
	$chr=$col[0];
	$str=$col[1];
	print OUT "$chr\t$str\t$fpREG{$key}\n";
}
close (OUT);

