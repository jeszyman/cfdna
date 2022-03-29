#Demultiplex fastq.gz samples using barcode affixed to read id

my $position = 1; #0 if barcode is on left side of sequence tag; 1 if on right
my $bases = 4; #size of barcode
my $num_samples = 12; #number of multiplexed samples
my $fastq = $ARGV[0]; #input fastq.gz
my $outdir = $ARGV[1]; #output directory
my $sample2bc = "";
if(@ARGV>2) {$sample2bc = $ARGV[2];} #demultiplex barcodes in sample2barcode.txt (optional)

#store barcodes
my %final_barcodes = ();

if(@ARGV == 2){ #sample2barcode.txt not specified

open(FILE, "zcat $fastq |") or die $!;

    my $itor = -1;

    my %barcodes = ();

    while(<FILE>){
        $itor++;
        #last if $itor >= 1000;
        my $line = $_;
        chomp($line);
        if($itor % 4 == 0){
        my @vars = split(":", $line);
        my $seq = $vars[@vars-1];
        my $index = 0;
        if($position == 1) {$position = length($seq)-$bases;}
        my $barcode = substr($seq,$position,$bases);
        #print "$seq\t$barcode\n";
        $barcodes{$barcode}++;
        }   
    }

    close(FILE);

    $itor = 0;
    foreach my $bc (sort {$barcodes{$b} <=> $barcodes{$a}} (keys %barcodes)){
        my $count = $barcodes{$bc};
        #print "$bc\t$count\n";
       my $out = $fastq;
        $out =~ s/.fastq.gz/_$bc.fastq/g;
        $out =~ s/.+\///g;
        $out = "$outdir/$out";
        open my $bcout,'>', $out or die "$!";
        $final_barcodes{$bc} = $bcout;
        #print $bc "$outdir/$out\n";
        $itor++;
        last if $itor == $num_samples;
    }

}else{
    
    open(FILE, $sample2bc) or die $!;
    while(<FILE>){
        my $line = $_;
        chomp($line);
        my @vars = split("\t", $line);
        my $seq = $vars[1];
        #if($position == 1) {$position = length($seq)-$bases;}
        #my $bc = substr($seq,$position,$bases);
	my $bc = $seq;
        $bc =~ s/N//g;
        my $out = $fastq;
        $out =~ s/.fastq.gz/_$bc.fastq/g;
        $out =~ s/.+\///g;
        $out = "$outdir/$out";
        open my $bcout,'>', $out or die "$!";
        $final_barcodes{$bc} = $bcout;
    }
}

#Final pass to print sequences by barcodes

open(FILE, "zcat $fastq |") or die $!;

my $barcode = "";

$itor = -1;

while(<FILE>){
    $itor++;
    #last if $itor >= 1000;
    my $line = $_;
    chomp($line);
    if($itor % 4 == 0){
        my @vars = split(":", $line);
        my $seq = $vars[@vars-1];
        my $index = 0;
        if($position == 1) {$position = length($seq)-$bases;}
        $barcode = substr($seq,$position,$bases);
        my @vars2 = split(" ", $line);
        if(@vars2 > 1) {$line = $vars2[0] . ":" . $seq . " " . $vars2[1];}
        }
        next if !(exists($final_barcodes{$barcode}));
        #print "$barcode\n";
        my $bcout = $final_barcodes{$barcode};
        print $bcout "$line\n";
}

#close output files
foreach my $bc (keys %final_barcodes){
    close($bc);
}
