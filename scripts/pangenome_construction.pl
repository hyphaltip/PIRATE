#!/usr/bin/env perl

# Pangenome construction using cd-hit to deflate clusters, all-vs-all BLAST, nested/single MCL clustering and reinflation.

use strict;
use warnings;

use Getopt::Long qw(GetOptions);
use Pod::Usage;
use File::Basename;
use File::Copy;
use FindBin;
use Cwd 'abs_path';
my $script_path = abs_path($FindBin::RealBin);

# Version

=head1  SYNOPSIS

	pangenome_construction.pl -i /path/to/fasta -o /path/to/output/

	Input/Output options:
	-i|--input	input fasta file [required]
	-o|--output	output directory [default: input directory]
	-l|--loci	file containing loci and genome as seperate columns 
			[required for core extraction during cdhit]
	
	Clustering options:
	-p|--perc	single % id threshold to use for pangenome 
			construction [default: 98]
	-s|--steps	multiple % id thresholds to use for pangenome 
			construction, comma seperated 
			[default: 50,60,70,80,90,95,98]
	-n|--nucl	create pangenome on nucleotide sequence 
			[default: amino acid]
	
	CDHIT options: 
	--cd-low	cdhit lowest percentage id [default: 98]
	--cd-step	cdhit step size [default: 0.5]
	--cd-core-off	don't extract core families during cdhit clustering 
			[default: on]
	--cd-mem	specify amount of memory required for CD-HIT in MB 
			[default: 5*input file size]
	
	BLAST options:
	-e|--evalue	e-value used for blast hit filtering [default: 1E-6]
	--diamond	use diamond instead of BLAST - incompatible 
			with --nucleotide [default: off]
	--diamond-split	split diamond files into batches for processing 
			[default: off] 	
	--hsp-prop	remove BLAST hsps that are < hsp_prop proportion
			of query length/query hsp length [default: 0]
	--hsp-len	remove BLAST hsps that are < hsp_len proportion
			of query length [default: 0]
	
	MCL options:
	-f|--flat	MCL inflation value [default: 1.5]
	
	General options:
	-r|--retain	do not delete temp files
	-t|--threads	number of threads/cores used to use [default: 2]
	-q|--quiet	switch off verbose
	-h|--help 	usage information

=cut

# switch off buffering
$| = 1; # turn off buffering for real time feedback.

# command line options
my $help = 0;
my $quiet = 0;
my $retain = 0;
my $nucleotide = 0;
my $threads = 2; 

my $input_file = '';
my $output_dir = '';
my $loci_list = '';

my $perc = 98;
my $steps = "50,60,70,80,90,95,98";
my $cd_low = 98;
my $cd_step = 0.5;
my $evalue = ""; # 0.001 for diamond
my $inflation_value = 1.5;
my $hsp_prop_length = 0.0;
my $hsp_length = 0;

my $cdhit_aS = 0.9;
my $m_required = 0;

my $diamond = 0;
my $diamond_split = 0;

my $exit_status = 1; 
my $core_off = 0;
my $cdhit_core = 1;

GetOptions(

	'help|?' 	=> \$help,
	'quiet'	=> \$quiet,
	'retain' => \$retain,
	
	'input=s' 	=> \$input_file,
	'output=s'	=> \$output_dir,
	'loci=s' 		=> \$loci_list,
	
	'threads=i'	=> \$threads,
	'steps=s'	=> \$steps,
	'perc=s'	=> \$perc,
	
	'cd-low=i' => \$cd_low,
	'cd-step=f' => \$cd_step,
	'cd-core-off' => \$core_off,
	'cd-hit-aS' => \$cdhit_aS,
	'cd-mem=i' => \$m_required,
	
	'flat=f' 	=> \$inflation_value,
	'evalue=f' => \$evalue,
	'hsp-prop=f' => \$hsp_prop_length,
	'hsp-length=f' => \$hsp_length, 
	
	'diamond' => \$diamond,
	'diamond-split' => \$diamond_split,

	'nucleotide' => \$nucleotide,
	
) or pod2usage(1);
pod2usage(1) if $help;
pod2usage if $input_file eq '';

# check dependencies - no version check.
my $cd_hit_bin = "";
my $cd_hit_est_bin = "";

my $dep_err = 0;
my $diamond_err = 0;
my $cdhit_check = 0;

# two alternative cd-hit invocations
if( (`command -v cdhit;`) && (`command -v cdhit-est;`) ){ 
	$cd_hit_bin = "cdhit";
	$cd_hit_est_bin = "cdhit-est";
	$cdhit_check = 1; 
}
if( (`command -v cd-hit;`) && (`command -v cd-hit-est;`) ){ 
	$cd_hit_bin = "cd-hit";
	$cd_hit_est_bin = "cd-hit-est";
	$cdhit_check = 1; 
}
print "cd-hit binary not found in system path.\n" if $cdhit_check == 0;
$dep_err = 1 if $cdhit_check == 0;

unless( `command -v blastp;` ){ 
	$dep_err = 1;
	print " - ERROR: blastp binary not found in system path.\n";
}
unless( `command -v blastn;` ){ 
	$dep_err = 1;
	print " - ERROR: blastn binary not found in system path.\n";
}
unless( `command -v makeblastdb;` ){ 
	$dep_err = 1;
	print " - ERROR: makeblastdb binary not found in system path.\n";
}
unless( `command -v mcl;` ){ 
	$dep_err = 1;
	print " - ERROR: mcl binary not found in system path.\n";
}
unless( `command -v mcxdeblast;` ){ 
	$dep_err = 1;
	print " - ERROR: mcxdeblast binary not found in system path.\n";
}
unless( (`command -v diamond makedb;`) && (`command -v diamond blastp;`) ){ 
	print " - WARNING: diamond binaries not found in system path.\n";
	$diamond_err = 1;
}
die " - ERROR: dependecies missing.\n" if $dep_err == 1;

# set cdhit core options
$cdhit_core = 0 if $core_off == 1;

# expand input and output files/directories
$input_file = abs_path($input_file);
my $input_dir = dirname(abs_path($input_file));
$output_dir = $input_dir if $output_dir eq '';

# make output directory if it doesn't exist. 
unless( -d "$output_dir" ){
	 die "could not make working directory in $output_dir\n" unless mkdir "$output_dir"; 
}
$output_dir = abs_path($output_dir);

# check for mandatory input arguments
pod2usage( {-message => q{input directory is a required arguement}, -exitval => 1, -verbose => 1 } ) if $input_dir eq ''; 

# check for presence of input/output directories.
pod2usage( {-message => "input directory:$input_dir is not a directory", -exitval => 1, -verbose => 1 } ) unless -d $input_dir; 
pod2usage( {-message => "ERROR: diamond can only be applied to protein alignments", -exitval => 1, -verbose => 1 } ) if ( ($diamond == 1) && ($nucleotide == 1) );
pod2usage( {-message => "ERROR: diamond binaries not found", -exitval => 1, -verbose => 1 } ) if ( ($diamond == 1) && ($diamond_err == 1) );

# working variables
my @files = ($input_file);
my $no_files = scalar(@files);
my $sample = "";

# identify thresholds and check they are numeric
my @thresholds = ();
if( $steps eq "" ){
	@thresholds = ($perc);
}
else{ 
	@thresholds = split (/\,/, $steps);
	@thresholds = sort  {$a <=> $b} @thresholds;
}
for (@thresholds) { if ( $_ !~ /^\d+$/ ) { die "Threshold value $_ is not numeric.\n" } }

# check prop is not > 1
die " - ERROR: hsp-prop is a proportion and must be < 1.\n" if $hsp_prop_length > 1; 
die " - ERROR: hsp-length is a proportion and must be < 1.\n" if $hsp_length > 1; 

# check files exist and have correct suffix.
my @suffix_list = (".aa.fasta" , ".fasta" , ".fa" , ".fas");
for my $file( @files ){

	my $suffix_check = 0;
	for (@suffix_list){ $suffix_check = 1 if $file =~ /$_$/  }
	die "$file suffix not recognised\n" if $suffix_check == 0;

	die "file $file does not exist.\n" if !-f $file;
		
} 

# check for conflicts between cd-hit low and %id threshold
my $max_thresh = $thresholds[scalar(@thresholds)-1];
die " - WARNING: lowest cd-hit threshold ($cd_low) is below blast % id value ($max_thresh)" if $cd_low < $max_thresh;

# parse and check loci list if passed via -l
my %loci;
my $no_loci = 0;
my $no_genomes = 0;
if( $loci_list ne '' ){
	$loci_list = abs_path($loci_list);
	die "loci_list file not found.\n" unless -f $loci_list;
	open LOCI, "$loci_list" or die $!;
	while (<LOCI>){	
		if(/^(\S+)\t(\S+)\t/){ 
			$loci{$1}=$2;
		}
	}	
	$no_loci = scalar ( keys(%loci) );
	my %no_g = map {$_ => 1} values(%loci);
	$no_genomes =  scalar( keys (%no_g) );
}

# set BLAST/DIAMOND e-value filter 
if ($evalue eq ""){
	if ( $diamond == 1 ){
		$evalue = "0.001";
	}else{
		$evalue = "1E-6";
	}
}

# user feedback
if ($quiet == 0 ){
	print "\nOptions:\n\n";
	print " - WARNING: cannot extract core loci during cdhit clustering unless loci list is provided!\n" if ( ($cdhit_core == 1) && ($loci_list eq '') );

	print " - Creating pangenome on nucleotide % identity.\n" if $nucleotide == 1;
	print " - Creating pangenome on amino acid % identity using BLAST.\n" if ( ($nucleotide == 0) && ($diamond == 0) );
	print " - Creating pangenome on amino acid % identity using DIAMOND.\n" if ( ($nucleotide == 0) && ($diamond == 1) );
	print " - Input directory:\t$input_dir\n";
	print " - Output directory:\t$output_dir\n";
	print " - Number of input files: $no_files\n";
	print " - Threshold(s): @thresholds\n";
	print " - MCL inflation value: $inflation_value\n";
	print " - Homology test cutoff: $evalue\n";
	print " - Recipocal length cutoff: $hsp_prop_length\n" if $hsp_prop_length > 0;
	print " - Loci file contains $no_loci loci from $no_genomes genomes.\n" if $loci_list ne '';
	print " - Extracting core loci during cdhit clustering\n" if ( ($cdhit_core == 1) && ($loci_list ne '') );
}

# make mcl temp dir
mkdir( "$output_dir/mcl_sub/" );

# variables for user feedback
my $processed = 0;
my $no_processed = 0;

# timer sub
my $time = time();
sub time_update {
	my $time_diff = time() - $time;
	print " - completed in $time_diff secs\n";
	$time = time();
}

# process files sequentially
for my $file( @files ){

	++$processed;	
	
	# find sample name
	my ($sample, $s_path, $suffix) = fileparse($file , @suffix_list);
				
	# feedback
	print " - Opening $sample\n" if $quiet == 0;
	
	# make temp fasta sequence without alignment characters and single line.
	my @seq = ();
	my @seq_out = ();
	my $header = "";
	my @g_lengths = ();
	my $temp_seq = "";
	my $seq_count = 0;
	my @seq_ids = ();
	my $no_included = "";
	
	open FILE, "$file" or die "$file not found.\n";
	while(<FILE>){
		
		my $line = $_;
		$line =~ s/\R//g;
				
		if(/^>/){
		
			#$seq_count++; 
						
			unless( scalar(@seq) == 0 ){
				$temp_seq = join( "" , @seq );
				push ( @seq_out, join("\n" , $header , $temp_seq) );
				push ( @g_lengths, length($temp_seq) );
				$seq_count++;
			}	

			$header = $line;	
			@seq = ();
			
			push(@seq_ids, $header); 
			
		}else{
		
			$line =~ s/\-//g;
			push( @seq, $line );
			
		}
	}close FILE;
	print " - $file contains $seq_count sequences.\n";
	
	# store final sample
	unless( scalar(@seq) == 0 ){
		$temp_seq = join( "" , @seq );
		push ( @seq_out, join("\n" , $header , $temp_seq) );
		push ( @g_lengths, length($temp_seq) );
		$seq_count++;
	}	
	
	# sort file on length.
	my @idx = sort { $g_lengths[$b] <=> $g_lengths[$a] } 0 .. $#g_lengths;
	
	# Print to all sequences file.
	open TEMP, ">$output_dir/$sample.all_sequences.fasta" or die $!;
	print TEMP join("\n", @seq_out[@idx]);
	close TEMP;
	
	# sanity checks - number of sequences in file.
	my $no_sequences = scalar ( @seq_out );
	die " - ERROR: No sequences in input fasta file.\n" if $seq_count == 0;
	die " - ERROR: Number of sequences ($no_sequences) do not match number of headers ($seq_count)\n" if $seq_count != $no_sequences;
	
	# cluster variables using CD-Hit.
	my %cluster_hash = ();
	
	# make cd-hit log file.
	my $cdhit_log = "$output_dir/$sample.cdhit_log.txt";
	open CD_LOG, ">$cdhit_log" or die $!;
	
	# create core file and core hash.
	my %core = ();
	open CORE, ">$output_dir/$sample.core_clusters.tab" or die "couldn't open $sample.core_clusters.tab\n";
	
	# make temporary fasta file for cd-hit
	if ( -f "$output_dir/$sample.all_sequences.fasta" ){
	 	copy( "$output_dir/$sample.all_sequences.fasta" , "$output_dir/$sample.temp.fasta") or die " - ERROR: could not create $output_dir/$sample.temp.fasta" ;
	}else{
		die " - ERROR: $output_dir/$sample.all_sequences.fasta was not found - did you run out of memory?" ;
	}
	
	my $no_reduced = 0;
	
	# run cd-hit at multiple thresholds.
	my $final_threshold = "";
	for (my $i = 100; $i >= $cd_low; $i -= $cd_step) {	   		 
	
		my $curr_thresh = $i/100;
		
		# Number of loci to pass to cd-hit
		$no_included = $no_sequences - (keys %core);
		
		# Filter core loci from temp fasta file.
		if ( keys(%core) > 0 ){	
		
			my $header = "";
			my $sample_check = 0;	
			
			open FASTA_IN, "$output_dir/$sample.temp.fasta" or die "$file not found.\n";
			open FASTA_OUT, ">$output_dir/$sample.temp2.fasta" or die "$file not found.\n";
			
			while(<FASTA_IN>){
		
				my $line = $_;
				$line =~ s/\R//g;
				
				if(/^\>(\S+)*/){
			
					$header = $1;
				
				}else{
		
					if( !$core{ $header } ){
						print FASTA_OUT ">$header\n$line\n";
						++$sample_check;
					}
				}
			}
			
			close FASTA_IN;
			close FASTA_OUT;
		
			# make filtered file the working file
			move( "$output_dir/$sample.temp2.fasta" , "$output_dir/$sample.temp.fasta") or die " - ERROR: could not rename $output_dir/$sample.temp.fasta";
		
			# Sanity check
			die "Number of samples in $output_dir/$sample.temp.fasta ($no_included) does not match number of included loci ($sample_check).\n" if $sample_check != $no_included;
		
		}
		
		# calculate memory for cdhit if not specified
		if ($m_required == 0){
			$m_required = -s "$output_dir/$sample.temp.fasta"; # bytes
			$m_required *= 5; # triple
			$m_required = int($m_required/1000000); #MB
			$m_required = 2000 if($m_required < 2000); # set minimum
		}
		
		# run cdhit
		print " - Passing $no_included loci to cd-hit at $i%  \n" if $quiet == 0;
		my $n = "";
		#my $cd_hit_out = "";
		if( $nucleotide == 0 ){
		
			# select appropriate word size
			if ( $curr_thresh > 0.7 ){
				$n = 5;
			}elsif ( $curr_thresh > 0.6 ){
				$n = 4;
			}elsif ( $curr_thresh > 0.5 ){
				$n = 3;
			}elsif ( $curr_thresh > 0.4 ){
				$n = 2;
			}else{
				$curr_thresh = 0.4;
				$n = 2;
				print " - WARNING: cluster threshold ($curr_thresh) below recommended setting.";
				print " - Setting cluster threshold (-c) to 0.4 and word size (-n) to 2";
			}

			# run cd-hit
			my $cd_hit_command = "$cd_hit_bin -i $output_dir/$sample.temp.fasta -o $output_dir/$sample.$i -aS $cdhit_aS -c $curr_thresh -T $threads -g 1 -n $n -M $m_required -d 256 >> $cdhit_log";
			print " - command: \"$cd_hit_command\"\n";
			`$cd_hit_command`;
			

		}else{
		
			# select appropriate word size
			my $n = "";
			if ( $curr_thresh > 0.90 ){
				$n = 11;
			}elsif ( $curr_thresh > 0.90 ){
				$n = 9;
			}elsif ( $curr_thresh > 0.88 ){
				$n = 7;
			}elsif ( $curr_thresh > 0.85 ){
				$n = 6;
			}elsif ( $curr_thresh > 0.80 ){
				$n = 5;
			}else{
				$curr_thresh = 0.8;
				$n = 4;
				print " - WARNING: cluster threshold ($curr_thresh) below recommended setting.";
				print " - Setting cluster threshold (-c) to 0.8 and word size (-n) to 4";
			}
		
			# run cdhit est
			my $cd_hit_command = "$cd_hit_est_bin -i $output_dir/$sample.temp.fasta -o $output_dir/$sample.$i -aS $cdhit_aS -c $curr_thresh -T $threads -g 1 -n $n -M $m_required -d 256 -r 0 >> $cdhit_log";
			print " - command: \"$cd_hit_command\"\n";
			`$cd_hit_command`;
		}
		die " - ERROR: cdhit failed - the errors are logged at $cdhit_log\n" if $?;
		
		# variables
		my $c_name = "";
		%cluster_hash = ();
		
		# Store cluster loci
		open CLUSTER, "$output_dir/$sample.$i.clstr" or die $!;
		while(<CLUSTER>){
	
			# Add clustered loci to storage hash.
			if(/^>Cluster\s(\d+)*/){

				$c_name = $1 + 1;

			}elsif( /^\d+\s+(\d+)(aa|nt)\,\s+>(.+)\.\.\./ ){ 
			
				# sanity check 
				die "cdhit header error\n" if $c_name eq "";
				
				# cluster hash
				$cluster_hash {$3} = $c_name;
				
			}else{
				die "$_ did not match cd-hit format.\n";
			}
			
		}close CLUSTER;
		
		# Identify core loci (dosage of one per genome).
		if( ( $loci_list ne '' ) && ( $cdhit_core == 1 ) ){
		
			#print " - Extracting core loci...\n" unless $quiet == 1;
			
			# variables
			my %cluster_genomes = ();
			my %cluster_count = ();
			
			# count no of genes and no of genomes
			for my $l( keys %cluster_hash ){
		
				die "no genome found for loci $l\n" if !$loci{$l};
				$cluster_genomes { $cluster_hash{$l} }{ $loci{$l} }++;
				$cluster_count { $cluster_hash{$l} } { $l } = 1 ;
				
			}
			
			# Identify core clusters.
			for my $lc ( keys %cluster_genomes ){
				
				my $c_count = scalar(keys(%{$cluster_genomes{$lc}}));
				my $g_count = scalar(keys(%{$cluster_count{$lc}}));
				
				if( ( $c_count == $no_genomes ) && ( $g_count == $no_genomes) ){
					
					# add to core_clusters.tab file.					
					my $core_out = join ( "\t" , keys(%{$cluster_count{$lc}}));
					print CORE "$core_out\n";
					
					# add to core_hash.
					foreach ( keys %{$cluster_count{$lc}} ){
						$core{$_} = $lc;
					}	
				}
						
			}		
			
		}
		
		# store threshold
		$final_threshold = $i;
		
	}
	close CORE;	
	close CD_LOG;
		
	# Remove core loci from cluster hash.
	my %cluster_names = ();
	for my $l ( keys %cluster_hash ){
		$cluster_names{ $cluster_hash{$l} }{ $l } = 1 if !$core{$l};
		delete $cluster_hash{$l} if $core{$l};
	}
	
	# Make filtered representative fasta file.
	$no_included = $no_sequences - (keys %core);
	$header = "";
	my $final_check = 0;	
	my %rep = ();
			
	open FASTA_IN, "$output_dir/$sample.$final_threshold" or die "$file not found.\n";
	open FASTA_OUT, ">$output_dir/$sample.representative.fasta" or die "$file not found.\n";
	
	while(<FASTA_IN>){

		my $line = $_;
		$line =~ s/\R//g;
		
		if(/^\>(\S+)*/){
	
			$header = $1;
		
		}else{

			if( !$core{ $header } ){
			
				# print to file
				print FASTA_OUT ">$header\n$line\n";
				
				# add to representative hash 
				$rep{$header} = 1;
				
				++$final_check;
			}
		}
	}
	
	close FASTA_IN;
	close FASTA_OUT;

	# Sanity check
	my $n_clusters = scalar(keys %cluster_names);
	die "Number of samples in $output_dir/$sample.representative.fasta ($final_check) does not match number of representative loci( $n_clusters ).\n" if $final_check !=  $n_clusters;

	# Store clusters by cluster_name
	open REP_CLUSTER, ">$output_dir/$sample.cdhit_clusters" or die $!;
	
	my @cluster_line = ();
	for my $cn ( keys %cluster_names ){
		
		# store all loci in cluster
		@cluster_line = "$cn\t-";
		for my $cl (keys %{$cluster_names{$cn}} ){
			push( @cluster_line , $cl);
		}
		
		# print to file
		print REP_CLUSTER join("\t", @cluster_line), "\n";

	}close REP_CLUSTER;

	# Sanity check - fasta file sequences == number representative sequences
	my $representative_check = keys ( %cluster_names ); 
	my $fasta_check = `grep ">" < $output_dir/$sample.representative.fasta | wc -l`;
	die "Number of samples in $output_dir/$sample.representative.fasta ($fasta_check) does not match number of loci to include ($representative_check).\n" if $fasta_check != $representative_check;
	time_update();
	
	# Sanity check - fasta no loci == number of starting sequences.
	my $final_loci_no = scalar ( keys ( %cluster_hash ) );
	my $core_cluster_no = scalar ( keys %core );
	my $total_clusters = $core_cluster_no + $final_loci_no;
	die "number of sequences clustered via cd-hit ($total_clusters) does not match number of input sequences ($no_sequences)" if $total_clusters != $no_sequences;
	
	# User feedback
	print "\n - $core_cluster_no core loci (",(($core_cluster_no/$total_clusters)*100), "%)\n" if ( ($quiet == 0) && ($cdhit_core == 1) );
	print " - $final_loci_no non-core loci (",(($final_loci_no/$total_clusters)*100), "%)\n" unless $quiet == 1;
	
	# clear variables 
	%core = ();
	
	# user feedback
	print "\n - ", scalar(keys(%cluster_names)), " representative loci passed to blast.\n" if $quiet == 0;
	
	# all vs all blast
	my $blast_in = "$output_dir/$sample.representative.fasta";
	my $blast_out = "$output_dir/$sample.blast.output";
	my $mask_file = "$output_dir/$sample.mask";
	
	# use appropriate program
	if( $nucleotide == 0 ){

		if ( $diamond == 1 ){
		
			print "\n - running all-vs-all DIAMOND on $sample\n" if $quiet == 0;
			`diamond makedb --in $blast_in --db $output_dir/$sample.diamond_db 2>/dev/null`;
			
			# output format
			my $outfmt = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore";
			$outfmt = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" if ( ($hsp_prop_length > 0) || ($hsp_length > 0) );
			
			# set query coverage as %
			my $q_cov = $hsp_length*100;
			
			if ($diamond_split == 1){
			
				# split file and run in parallel 
				`cat $blast_in | parallel --recstart '>' --jobs $threads --pipe "diamond blastp -d $output_dir/$sample.diamond_db -c 1 --query-cover $q_cov --masking 0 --evalue $evalue --max-hsps 1 --threads 1 --outfmt $outfmt --sensitive --max-target-seqs 0" > $blast_out 2>/dev/null`;
				
			}else{

				# run as one file (faster than parallel)
				`diamond blastp -q $blast_in -d $output_dir/$sample.diamond_db -c 1 --query-cover $q_cov --masking 0 --evalue $evalue --max-hsps 1 --threads $threads --outfmt $outfmt --sensitive --max-target-seqs 0 > $blast_out 2>/dev/null`;
				
			}
			
			# error check
			die "diamond blastp failed.\n" if $?;
			
		}
		else{
		
			print "\n - running all-vs-all BLASTP on $sample\n" if $quiet == 0;
			
			# mask sequences using segmasker
			#`segmasker -in $blast_in -infmt fasta -outfmt maskinfo_asn1_bin -out $mask_file`;	
			#`makeblastdb -in $blast_in -dbtype prot -mask_data $mask_file`;
			
			# make db
			`makeblastdb -in $blast_in -dbtype prot`;
			
			# output format
			my $outfmt = "\'6\'";
			$outfmt = "\'6 std qlen slen\'" if ( ($hsp_prop_length > 0) || ($hsp_length > 0) );
			
			# run blastp
			`cat $blast_in | parallel --recstart '>' --jobs $threads --pipe "blastp -max_hsps 1 -outfmt $outfmt -num_threads 1 -max_target_seqs 10000 -evalue $evalue -query - -db $blast_in" > $blast_out`;
			die "blastp failed.\n" if $?;
			
		}
	
	}else{
	
		# output format
		my $outfmt = "\'6\'";
		$outfmt = "\'6 std qlen slen\'" if ( ($hsp_prop_length > 0) || ($hsp_length > 0) );
	
		print "\n - running all-vs-all BLASTN on $sample\n" if $quiet == 0;
		`makeblastdb -in $blast_in -dbtype nucl`;
		`cat $blast_in | parallel --recstart '>' --jobs $threads --pipe "blastn -task \"blastn\" -outfmt $outfmt -num_threads 1 -dust no -evalue $evalue -max_target_seqs 10000 -max_hsps 1 -query - -db $blast_in" > $blast_out`;
		die "blastn failed.\n" if $?; 
		
	}
	time_update();
	
	# ensure blast output has all representative sequences vs themselves (short sequences maybe removed on e-value).
	open BLAST_OUT, "$blast_out" or die $!;
	open BLAST_TEMP, ">$blast_out.temp" or die $!;
	#my $av_bit = ""; # average bit score of same-same blast results for dataset.
	my $max_bit = ""; # minimum bit score of same-same blast results for dataset.
	while (<BLAST_OUT>){
	
		my $line = $_;
		my @line = split(/\t/, $line, -1);
	
		# identify same-same lines
		if( $line[0] eq $line[1] ){
		
			# same-same comparison exists - do not print 
			#$rep { $line[0] } = 1; # print placeholder line
			$rep { $line[0] } = 2; # do not print placeholder line
			
			# find max bit score
			$max_bit = $line[11] if $max_bit eq "";
			$max_bit = $line[11] if $line[11] < $max_bit;
			
			# ensure comparison has 100% identity
			$line[2] = "100.0";			
			
			# print line 
			print BLAST_TEMP join("\t", @line);
			
			# removed - used to use average and print existing lines
			#$av_bit = $line[11] if $av_bit eq "";
			#$av_bit = ($av_bit + $line[11]) / 2;
		
		}
		# [optional] filter on hsp length as a proportion of the total length of the quesry sequence
		elsif( $hsp_length > 0 ){
	
			my $q_len = $line[12];
			my $q_hsp_len = ($line[7] - $line[6]) + 1;
		
			# test hsp length against original sequence length and subject hsp alignment length vs query hsp length.			
			my $hspVSq = $q_hsp_len / $q_len;
			$hspVSq = 1 - ($hspVSq - 1) if $hspVSq > 1; 
		
			# print if both are > hsp_prop_length
			if ( ( $hspVSq > $hsp_length ) ){
				print BLAST_TEMP "$line";
			}
				
		}
		# [optional] filter on hsp length hsp percentage length < hsp_prop_length removed.
		elsif( $hsp_prop_length > 0 ){
		
			my $q_len = $line[12];
			my $q_hsp_len = ($line[7] - $line[6]) + 1;
			my $s_hsp_len = ($line[9] - $line[8]) + 1;
			
			# test hsp length against original sequence length and subject hsp alignment length vs query hsp length.			
			my $hspVSq = $q_len / $s_hsp_len;
			$hspVSq = 1 - ($hspVSq - 1) if $hspVSq > 1; 
			my $hspVShsp = $q_hsp_len / $s_hsp_len;
			$hspVShsp = 1 - ($hspVShsp - 1) if $hspVShsp > 1; 
			
			# print if both are > hsp_prop_length
			if ( ( $hspVSq > $hsp_prop_length ) && ( $hspVShsp > $hsp_prop_length) ){
				print BLAST_TEMP "$line";
			}
					
		}
		# otherwise print
		else{
		
			print BLAST_TEMP "$line";
			
		}
		
	}
	close BLAST_OUT;
	close BLAST_TEMP;
	
	# replace original file.
	move( "$blast_out.temp" , "$blast_out") or die " - ERROR: could not rename $blast_out.temp";
	
	# Add same-same hits for each representative sequence - this ensures none are lost by filtering/MCL.
	open BLAST_OUT, ">>$blast_out" or die $!;
	for my $rloci ( keys %rep ){

		if ( $rep{$rloci} == 1){
		
			# uses a placeholder value for a same-same sequence - previously used av_bit/max bit
			
			if ( ($hsp_prop_length > 0) || ($hsp_length > 0) ){
				print BLAST_OUT "$rloci	$rloci	100	1234	0	0	1	1234	1	1234	0	2335	1234	1234\n";
			}else{
				print BLAST_OUT "$rloci	$rloci	100	1234	0	0	1	1234	1	1234	0	2335\n";
			}
			#print BLAST_OUT "$rloci	$rloci	100.00	100	0	0	1	100	1	100	0.0e+00	$max_bit\n";
		} 
		
	}
	close BLAST_OUT;
		
	# clear representative sequence variable 
	%rep = ();
	
	# make mcl log file
	my $mcl_log = "$output_dir/$sample.mcl_log.txt";
	`echo -n "" > $mcl_log`;
	
	# Filter BLAST files at all thresholds and perform MCL on filtered hits.
	# Iterate through all clusters at higher thresholds.
	#print "\n - running mcl on $sample\n" if $quiet == 0;
	print "\n" if $quiet == 0;
	my $previous_clusters = ""; 
	for my $c(0..(scalar(@thresholds)-1) ){
	
		my $threshold = $thresholds[$c];
		
		# feedback
		print " - running mcl on $sample at $threshold    \n" if $quiet == 0;
		
		# filter blast results on threshold.
		`awk '{if(\$3 >= $threshold){print \$0}}' < $blast_out > $output_dir/$sample.$threshold.blast`;
				
		if( $threshold ==  $thresholds[0] ){
	
			# reformat to abc and run mcl on bitscores normalised by hsp length
			`mcxdeblast --line-mode=abc --m9 --score=r $output_dir/$sample.$threshold.blast 2>> $mcl_log | mcl - --abc -te $threads -I $inflation_value -o $output_dir/$sample.mcl_$threshold.clusters 2>> $mcl_log`;
			die "mcl failed at $threshold" if $?;
			
			# set working file for next iteration
			$previous_clusters = "$output_dir/$sample.mcl_$threshold.clusters";
			
			# feedback 
			my $no_clusters = `cat $output_dir/$sample.mcl_$threshold.clusters | wc -l`;
			$no_clusters =~ s/\n//;
			print " - $no_clusters clusters at $threshold %";
		
		}else{
		
			# make mcl cluster file.
			`echo -n "" > $output_dir/$sample.mcl_$threshold.clusters`;
			
			# run mcl on each subcluster independently.
			my $cluster_no = 0;
			my %l_hash = ();
			my $ct_file = "$output_dir/$sample.mcl_$thresholds[$c-1].clusters";
			open CLUSTERS, "$ct_file" or die $!;
			while (<CLUSTERS>){
			
				$cluster_no++;
			
				my $line = $_;
				chomp $line;
								
				my @loci = split(/\s+/, $line);
				foreach( @loci ){ $l_hash{$_} = $cluster_no };
				
			}close CLUSTERS;
			
			# check for clusters from previous iteration
			die " - ERROR: no clusters in $ct_file\n" if $cluster_no == 0;
			
			# make empty file for filtered blast results.
			foreach(1..$cluster_no){ 
				open BLAST_SUB, ">$output_dir/mcl_sub/cluster_$_.blast" and close BLAST_SUB or die $!;
			}
			
			# filter blast output into multiple temp cluster files.
			open BLAST, "$output_dir/$sample.$threshold.blast" or die $!;
			while(<BLAST>){
			
				my $b_line = $_;
			
				if(/^(\S+)\s+(\S+)\s+/){					
								
					if( $l_hash{$1} == $l_hash{$2} ){
						open BLAST_SUB, ">>$output_dir/mcl_sub/cluster_$l_hash{$1}.blast" or die $!;
						print BLAST_SUB $b_line ;
						close BLAST_SUB;
					}

				}
				
			}close BLAST;
			
			# make file containing blast file location and output file parallel mcl.
			open TEMP, ">$output_dir/mcl_sub/list.txt" or die $!;
			foreach(1..$cluster_no){ 
				print TEMP "$output_dir/mcl_sub/cluster_$_.blast\t$output_dir/mcl_sub/$sample.mcl_$_.clusters\n" or die $!;
			}close TEMP;			
			
			# run mcl in parallel. 
			`parallel -a $output_dir/mcl_sub/list.txt --jobs $threads --colsep '\t' \"mcxdeblast --line-mode=abc --m9 --score=r {1} 2>>$mcl_log | mcl - --abc -te 1 -I $inflation_value -o {2} 2>>$mcl_log \"`;
			die " - ERROR: mcl failed at $threshold\n" if $?;
			
			# compile clusters into one file for next iteration and remove original mcl cluster file.
			open MCL_CONCAT, ">$output_dir/$sample.mcl_$threshold.clusters" or die $!;
			for my $j (1..$cluster_no){ 
			
				open MCL_CLUSTER, "$output_dir/mcl_sub/$sample.mcl_$j.clusters" or die " - ERROR: no MCL cluster file for $sample - $j at $threshold\n";
				while (<MCL_CLUSTER>){
					print MCL_CONCAT $_;
				}close MCL_CLUSTER;
				
				# Remove temp mcl files.
				unlink "$output_dir/mcl_sub/cluster_$j.blast";
				unlink "$output_dir/mcl_sub/$sample.mcl_$j.clusters";
				
			}close MCL_CONCAT;
			
			# set working file for next iteration
			$previous_clusters = "$output_dir/$sample.mcl_$threshold.clusters";
			unlink "$output_dir/mcl_sub/list.txt";
			
			# feedback 
			my $no_clusters = `cat $output_dir/$sample.mcl_$threshold.clusters | wc -l`;
			$no_clusters =~ s/\n//;
			print " - $no_clusters clusters at $threshold %";
						
		}
		time_update();
	}
	print "\n" if $quiet == 0;
	
	# Reinflate clusters
	print " - reinflating clusters for $sample" if $quiet == 0;
	for my $t(@thresholds){
		
		# open output
		open INFLAT, ">$output_dir/$sample.$t.reclustered.reinflated" or die $!;
		
		# open input
		open CR, "$output_dir/$sample.mcl_$t.clusters" or die $!;

		# sanity check on number of output sequences.
		my $seq_count = 0;
		
		# loop through all clusters and reinflate where appropriate.
		while (<CR>){
			
			my $line = $_;
			$line =~ s/\R//;
			
			my @clusters_reinflated = ();

			foreach my $inflat( split(/\t/ , $line) ){
				
				# if cluster was previously deflated with cd-hit add in missing loci.				
				die "\nERROR: $inflat not in a cdhit cluster.\n" if !$cluster_hash{$inflat} ;
					
				my $index  = $cluster_hash{$inflat};	
				foreach( keys %{$cluster_names{$index}} ){
					push( @clusters_reinflated, $_);
					++$seq_count;
				}					
				
			}
			
			# print to file.
			print INFLAT join ("\t", @clusters_reinflated), "\n";			

		}				
		close INFLAT;
		close CR;
		
		# add core clusters
		`cat $output_dir/$sample.core_clusters.tab $output_dir/$sample.$t.reclustered.reinflated > $output_dir/temp.txt`;
		move("$output_dir/temp.txt" , "$output_dir/$sample.$t.reclustered.reinflated") or die " - ERROR: could not rename $output_dir/temp.txt";
		
		# check number of clusters in final file matches input.
		die "Reinflated sequences (", ($seq_count+$core_cluster_no) , ") does not match input number of sequences ($no_sequences) at $t threshold in sample $sample.\n" if ($seq_count+$core_cluster_no) != $no_sequences;
		
	}
	
	# clean up temporary files
	unless ( $retain == 1 ){
	
		rmdir "$output_dir/mcl_sub/";
		unlink "$output_dir/$sample.all_sequences.fasta";
		unlink "$output_dir/$sample.cdhit_log.txt";
		#unlink "$output_dir/$sample.blast.output";
		unlink "$output_dir/$sample.temp.fasta";
		
		if ( $nucleotide == 0 ){
		
			if ( $diamond == 1 ){
				unlink "$output_dir/$sample.diamond_db.dmnd";
			}else{
				unlink "$output_dir/$sample.representative.fasta.phr";
				unlink "$output_dir/$sample.representative.fasta.pin";
				unlink "$output_dir/$sample.representative.fasta.psq";
			}
			
		}else{
			unlink "$output_dir/$sample.blast.input.nhr";
			unlink "$output_dir/$sample.blast.input.nin";
			unlink "$output_dir/$sample.blast.input.nsq";
		}
		
		if ( $diamond == 1 ){
			unlink "$output_dir/$sample.blast.input.nsq";
		}
	
		for (@thresholds){ 
			unlink "$output_dir/$sample.$_.blast";
			unlink "$output_dir/$sample.$_.temp_sub.blast";
			unlink "$output_dir/$sample.mcl_$_.clusters";
			unlink "$output_dir/$sample.mcl_$_.temp.clusters" if -f "$output_dir/$sample.mcl_$_.temp.clusters";
		}
	
		for (my $j = 100; $j >= $cd_low; $j -= $cd_step) {	
			unlink "$output_dir/$sample.$j";
			unlink "$output_dir/$sample.$j.clstr";
	
		}
		
	}
	
}

# user feedback
print "\n - Finished\n\n" if $quiet == 0;

exit
