#script to process raw paired-end fastq.gz reads with 5' molecular barcodes

#=====================================================================
use Compress::Zlib;
use Cwd 'abs_path';
use File::Basename;

my $debug=0;
my $fastq = $ARGV[0]; #fastq(.gz) read 1 file
my $fastq2 = $ARGV[1]; #fastq(.gz) read 2 file
my $file;

# Will extract long/short barcodes based on demux ID unless forced:
# Force long == 1.
# Force short == 2.
my $force_length = $ARGV[2]; #Force long adapters if defined.
my $long=1;

my %demux_long = (
# Long demux IDs = 1
        AGGT => 1,
        CACA => 1,
        GCTA => 1,
        TGTG => 1,
        CTTC => 1,
        TCCT => 1,
# Short demux IDs = 0
        ATCG => 0,
        GTGT => 0,
        TAGC => 0,
        ACAC => 0,
        CGAT => 0,
        GAAG => 0
);

# Accepted punctuation marks
my %valid_puncts = (
# Long barcode puncts
        1 => ["GT", "CT"],
# Short barcode puncts
        0 => ["GT"]
);

my $rightwobble = 0; #wiggle room from the right of the fixed sequence for a match
my $barcodelen = 4; #length of inserted barcode (left side of sequence)

my %passed = (); #store sequences with legit barcodes from read 1 data

# Hashes for punctuation/barcode quality tracking
my %puncts = (); #track punctuation marks
my %puncts_q1 = ();
my %puncts_q2 = ();
my %barcodes = (); #track counts of each barcode
my %barcodes_q1 = ();
my %barcodes_q2 = ();
my %barcodes_nopunct = (); #Barcode region w/ invalid punct
my %barcodes_nopunct_q1 = ();
my %barcodes_nopunct_q2 = ();

#Extract in 3 iterations:
# i==0)  filter read 1 data and store passing sequence numbers and barcodes
# i==1)  print matching reads from read 2 data and combine barcodes in order (read 1 barcode: read 2 barcode)
# i==2)  print matching reads from read 1 data using combined barcodes in same order

for (my $i = 0; $i < 3; $i++) {

    if($i==1) {$file = $fastq2;}
    else {$file = $fastq;}

    if ($file =~ /\.gz/) { 
        open(FILE, "gunzip -c $file |") or die $!;
    } else { open(FILE, $file) or die $!; }

    my $linecounter = 0;
    my $header = "";
    my $punct = "";
    my $barcode = "";
    my $megabarcode="";
    my $printseq = 0;
    my $seqcounter = 0;
    my $outputseqitor = 0;

    my $output = $file;

    if ($file =~ /\.gz/) {
        $output =~ s/fastq\.gz/fastq/g;
    }
    $output =~ s/\_1\_/\_R1\_/g;
    $output =~ s/\_2\_/\_R2\_/g;

    if($i>0) {open(out, ">$output") or die $!;}

     while(<FILE>){
	my $line = $_;
	chomp($line);

        # FASTQ header line
        if($linecounter % 4 == 0){
            $header = $line;
            $header =~ s/ /:/g;
            # Choose long/short barcode based on demux ID - default long.
            if ($force_length == 1) {$long=1;$barcodelen=4;}
            elsif ($force_length ==2) {$long=0;$barcodelen=2;} 
            else {
                my $demux_id = substr($header, length($header)-4, 4);
                if ($demux_long{$demux_id}) {$long=1;$barcodelen=4;}
                else {$long=0;$barcodelen=2;}
            }
            # printseq==1 if punctuation mark is valid for matching R1/R2 reads. 
            $printseq = 0;
            $seqcounter++;
        }
        # FASTQ sequence line
        elsif($linecounter % 4 == 1){

            $punct=substr($line,$barcodelen,2);
            # Log punctuation mark frequency
            if ($i>0) { $puncts{$punct} += 1; }

            # Check punctuation mark validity
            my $punct_valid=0;
            foreach(@{$valid_puncts{$long}}) {
                if (substr($line, $barcodelen, length($_)-$rightwobble) eq substr($_, 0, length($_)-$rightwobble)) { $punct_valid=1; }
            }

            # Extract barcode, process if punctuation is valid.
            $barcode = substr($line,($barcodelen-2), 2);
            if ($punct_valid) {
                # First pass
                if($i == 0) { $passed{$seqcounter} = $barcode; }
                # Second/third passes
                elsif (exists($passed{$seqcounter})) {
                    $megabarcode = $passed{$seqcounter} . ":" . $barcode;
                    if($i==2) {$megabarcode = $barcode . ":" . $passed{$seqcounter};}
                    $passed{$seqcounter} = $barcode;
                    $barcodes{$barcode} += 1;
                    $outputseqitor++;
                    $printseq = 1;

                    #Exclude barcodes with N-bases in either R1/R2 
                    if (!($megabarcode =~ /N/)) {
                        print out  "$header:$megabarcode\n" . substr($line, $barcodelen + length($punct), length($line)) . "\n";
                    }
                }
             }

            if ($i>0 && !$printseq) {
                    $barcodes_nopunct{$barcode} += 1;
                     delete $passed{$seqcounter};
             }

        } elsif($linecounter % 4 == 2 & $i>0){
            # FASTQ "+" line
            if($printseq && !($megabarcode =~ /N/)) {print out  "$line\n";}
        } elsif ($linecounter % 4 == 3 & $i>0) {
            # FASTQ quality line
            # Quality scores for barcode/punctuation
            $puncts_q1{$punct} += unpack("C*",substr($line,$barcodelen,1))-33;
            $puncts_q2{$punct} += unpack("C*",substr($line,$barcodelen+1,1))-33;
            if ($printseq == 1) {
                $barcodes_q1{$barcode} += unpack("C*",substr($line,($barcodelen-2),1))-33;
                $barcodes_q2{$barcode} += unpack("C*",substr($line,($barcodelen-1),1))-33;    
                if (!($megabarcode =~ /N/)) {print out substr($line, $barcodelen + length($punct), length($line)) . "\n";}
            } else {
                $barcodes_nopunct_q1{$barcode} += unpack("C*",substr($line,($barcodelen-2),1))-33;
                $barcodes_nopunct_q2{$barcode} += unpack("C*",substr($line,($barcodelen-1),1))-33;
            }
        }
        $linecounter++;

    }

    if($i>0) {close(out);}
}

my $logname = $fastq;
$logname =~ s/\.fastq(\.gz)?/\.barcodeQC\.log/g;
$logname =~ s/\_1\_/\_/g;
print "Log file: $logname\n";
open(QCLOG, ">$logname") or die $!;

# Print freqs for barcode/punct combos
print QCLOG "Punctuation distribution and per-base mean quality:\n";
print QCLOG "Mark\t\%\tCount\tB1Q\tB2Q\n";
my $puncts_n=0; my $puncts_valid=0; my $puncts_total=0;
my @puncts_keys = sort { $puncts{$b} <=> $puncts{$a} } keys %puncts ;
foreach(@{$valid_puncts{$long}}) {
        $puncts_valid += $puncts{$_};
}
foreach(@puncts_keys) {
    $puncts_total += $puncts{$_};
}
foreach(@puncts_keys) {
    if ($_ =~ /N/) { $puncts_n += $puncts{$_}; }
    printf QCLOG "$_\t%.2f\t$puncts{$_}\t%.1f\t%.1f\n", ($puncts{$_}*100)/$puncts_total, $puncts_q1{$_}/$puncts{$_}, $puncts_q2{$_}/$puncts{$_};
}
print QCLOG "Valid punctuation: " . join(", ", @{$valid_puncts{$long}}) . "\n";
print QCLOG "Total reads: $puncts_total\n";
print QCLOG "Valid punctuation: $puncts_valid\n";
print QCLOG "N-containing punctuation: $puncts_n\n";
printf QCLOG "Percentage valid punctuation: %.2f\n", ($puncts_valid*100) / $puncts_total;
printf QCLOG "Percentage N-content punctuation: %.2f\n\n", ($puncts_n*100) / $puncts_total;

print QCLOG "Accepted barcode distribution and per-base mean quality:\n";
print QCLOG "Barcode\t\%\tCount\tB1Q\tB2Q\n";
my @barcodes_keys = sort { $barcodes{$b} <=> $barcodes{$a} } keys %barcodes ;
my $barcodes_n=0; my $barcodes_total=0;
foreach(@barcodes_keys) {
    $barcodes_total += $barcodes{$_};
}
foreach(@barcodes_keys) {
    if ($_ =~ /N/) { $barcodes_n += $barcodes{$_}; }
    printf QCLOG "$_\t%.2f\t$barcodes{$_}\t%.1f\t%.1f\n", ($barcodes{$_}*100)/$barcodes_total, $barcodes_q1{$_}/$barcodes{$_}, $barcodes_q2{$_}/$barcodes{$_};
}
print QCLOG "Total accepted barcodes: $barcodes_total\n";
print QCLOG "N-containing barcodes: $barcodes_n\n";
printf QCLOG "Percentage N-content accepted barcodes: %.2f\n\n", ($barcodes_n*100) / $barcodes_total;

print QCLOG "Unused barcode distribution and per-base mean quality:\n";
print QCLOG "Barcode\t\%\tCount\tB1Q\tB2Q\n";
my @barcodes_nopunct_keys = sort { $barcodes_nopunct{$b} <=> $barcodes_nopunct{$a} } keys %barcodes_nopunct ;
my $barcodes_nopunct_n=0; my $barcodes_nopunct_total=0;
foreach(@barcodes_nopunct_keys) {
    $barcodes_nopunct_total += $barcodes_nopunct{$_};
}
foreach(@barcodes_nopunct_keys) {
    if ($_ =~ /N/) { $barcodes_nopunct_n += $barcodes_nopunct{$_}; }
    printf QCLOG "$_\t%.2f\t$barcodes_nopunct{$_}\t%.1f\t%.1f\n", ($barcodes_nopunct{$_}*100)/$barcodes_nopunct_total, $barcodes_nopunct_q1{$_}/$barcodes_nopunct{$_}, $barcodes_nopunct_q2{$_}/$barcodes_nopunct{$_};
}
print QCLOG "Total unused barcodes: $barcodes_nopunct_total\n";
print QCLOG "N-containing unused barcodes: $barcodes_nopunct_n\n";
printf QCLOG "Percentage N-content unused barcodes: %.2f\n", ($barcodes_nopunct_n*100) / $barcodes_nopunct_total;
close(QCLOG);
