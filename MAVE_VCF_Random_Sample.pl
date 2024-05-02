# random sample n number of sites of a vcf based on a given percentage
# input vcf file or vcf.gz
# Author: Miguel Vallebueno  CC BY-NC-SA 4.0
# Date: 2024-05-2


##########example how to run:
# perl MAVE_VCF_Random_Sample.pl input_file.vcf percentage

$file=$ARGV[0];
$prp=$ARGV[1];

if($file =~ /$vcf.gz/){open(FIL,"gunzip -c $file |") or die; $file=~s/\.vcf.gz//;}
elsif($file =~ /$.vcf/){open(FIL,"$file"; $file=~s/\.vcf//;) or die;}
else{print "Error 1: your file does not end in .vcf or .vcf.gz";}


$out="$file\_smp$prp.vcf";


open(OUT,">$out") or die;

$cc=0;
while($line=<FIL>) {
    if($line =~ /^#/){print OUT "$line";next;}
    $cc++;
    my $rs = int(rand(100));  ### random sampler
    if ($rs > $prp) {next;}
    print OUT "$line";
}
close(OUT);
close(FIL);
print "TOTAL NUMBER of SITES <$cc>\n";
