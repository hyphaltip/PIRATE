#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long qw(GetOptions);
use Pod::Usage;
use File::Basename;
use Bio::Seq;
use Cwd 'abs_path';

# select representative sequences from each gene family

# Usage:

=head1  SYNOPSIS

 select_representative -i /path/to/PIRATE.gene_families.tab -g /path/to/gff/files/

 -i|--input		input PIRATE.gene_families.tab file [required]
 -g|--gff		gff file directory [required]
 -o|--output		output root, ffn and faa files will be printed [required]	
 -t|--threshold		% threshold below which clusters are excluded
 			[default: off]
 -m|--max-threshold	% threshold above which clusters are excluded
 			[default: off]
 -d|--dosage		upper threshold of dosage to exclude from alignment 
 			[default: off]
 -a|--alleles		store sequence for alleles not gene families 
 			[default: off]
 -q|--quiet		verbose off [default: on]
 -h|--help		usage information

=cut

# switch off buffering
$| = 1; # turn off buffering for real time feedback.

# command line options
my $help = 0;
my $quiet = 0;
my $aa = 0;
my $alleles = 0;

my $threshold = 0;
my $max_threshold = 100;
my $dosage_threshold = 0;

my $input_file = '';
my $gff_dir = '';
my $output_dir = '';

GetOptions(

	'help|?' 	=> \$help,
	'quiet' 	=> \$quiet,
	'input=s' 	=> \$input_file,
	'gff=s' 	=> \$gff_dir,
	'output=s'	=> \$output_dir,
	'threshold=f'	=> \$threshold,
	'max-threshold=f' => \$max_threshold,
	'dosage=f'	=> \$dosage_threshold,
	'alleles' => \$alleles,
	
) or pod2usage(1);
pod2usage(1) if $help;

# make paths absolute
$input_file = abs_path($input_file);
$gff_dir = abs_path($gff_dir);

# expand input and output files/directories
die "Input file not found.\n" unless -f "$input_file";
my $input_dir = dirname(abs_path($input_file));

# chack gff directory exists.
die " - ERROR: GFF directory not found.\n" unless -d "$gff_dir";

# check for mandatory input arguments
pod2usage( {-message => q{input directory is a required argument}, -exitval => 1, -verbose => 1 } ) if $input_dir eq ''; 
pod2usage( {-message => q{output directory is a required argument}, -exitval => 1, -verbose => 1 } ) if $output_dir eq ''; 

# Group/Loci variables
my %loci_group; # group of loci
my %group_list; # all loci in group
my %products;
my %gene_names;
my %group_names;

my @headers = ();
my @genomes = ();
my $total_genomes = 0;
my $no_headers = 0;
my $idx = 19; 

# Parse all groups and exclude gene families that do not meet threshold requirements
open GC, "$input_file" or die "$!";
while(<GC>){
	
	my $line =$_;
	$line =~ s/\R//g;
	
	if(/^allele_name\t/){
		
		# select index colum for input type/version
		$idx = 20 if $line =~ /\tno_loci\t/;
		$idx = 22 if $line =~ /\tcluster_order\t/;
		
		@headers = split (/\t/, $line, -1 );
		$no_headers = scalar(@headers);
		
		@genomes = @headers[$idx.. ($no_headers-1) ];
		$total_genomes = scalar(@genomes);		
		
	}else{
		
		# sanity check 
		die "No header line in $input_file" if scalar(@headers) == 0;
		
		my @l = split (/\t/, $line, -1 );
		
		# define group values
		my $allele = $l[0];
		my $group = $l[1];
		my $no_genomes = $l[6];
		my $dosage = $l[7]; # max dosage
		my $per_genomes = ($no_genomes / $total_genomes) * 100;
		my $prod = $l[3];
		my $gene = $l[2];
		
		# filter on thresholds
		if ( ($per_genomes >= $threshold) && ($per_genomes <= $max_threshold) && !($group_list{ $group }) ){
		
			# [optional] filter on dosage 
			if ( ($dosage_threshold == 0) || ($dosage <= $dosage_threshold) ){
			
				# store product and gene info
				if ($alleles == 1){
					$products {$allele} = $prod;
					$gene_names{$allele} = $gene;
				}else{
					$products {$group} = $prod;
					$gene_names{$group} = $gene;
				}
				
				# Store all loci for group
				for my $idx ( $idx..$#l ){
			
					my $entry = $l[$idx];
					my $entry_genome = $headers[ $idx ];
					
					unless( $entry eq "" ){
				
						$entry =~ s/\(|\)//g;
			
						foreach my $split_entry ( split(/;|:|\//, $entry, -1) ){
						
							# store alleles
							if ( $alleles == 1 ){
								$loci_group { $split_entry } {$allele} = 1;
								$group_list { $allele } = $no_genomes;
								$group_names { $allele } = $group;
							}
							else{
								$loci_group { $split_entry } {$group} = 1;
								$group_list { $group } = $no_genomes;
							}
							
						}								
				
					}
		
				}
				
			}
		
		}
	
	}
	
}

# Feedback
my $no_groups =  scalar ( keys %group_list );
print " - number of groups : $no_groups\n" if $quiet == 0;

# Check from presence gffs for all genomes.
for my $g ( @genomes ){
	die "No gff file for $g in $gff_dir\n" unless -f "$gff_dir/$g.gff";
}

# variables
my %stored_groups;
my %stored_length = ();
my %stored_seq = ();

# Extract sequence for all loci - only keep longest sequence for each group.
print " - extracting sequences from gffs\n" if $quiet == 0;
for my $genome ( @genomes ){
	
	# Extract sequences and print to file.
	my $count = 0;
	my $include = 0;
	my $contig_id = "";
	my @c_seq = ();
	my %seq_store;
	
	open INPUT, "$gff_dir/$genome.gff" or die $!;
	while(<INPUT>){
	
		my $line=$_;
		chomp $line;

		# start storing sequence after fasta header.
		if($line =~ /^##FASTA/){ 
			$include = 1; 
		}
		elsif( $include == 1){
		
			# header is used as contig id.
			if($line =~ /^>(.+)$/){
			
				++$count;
				
				# Don't store on first header
				if( $count > 1 ){
					my $store_seq = join ("", @c_seq);
					$seq_store{$contig_id} = $store_seq;
				}
				
				# Set variables
				$contig_id = $1;
				@c_seq = ();
						
			}
			
			# sequence is stored in array.
			else{ 
			
				# remove whitespace
				$line =~ s/\s+//g;
				
				if($line =~ /^([ATGCNatcgn-]+)$/){

					# sanity check - each contig should have id
					die "Contig has no header" if $contig_id eq "" ;
			
					# store uppercase sequence.
					my $seq = $1;
					$seq = uc($seq);
				
					push(@c_seq, $seq);
		
				}else{
				
					# replace characters that are not known bases with Ns.
					my $rep_chars = $line;
					$rep_chars =~ s/[ATGCNatcgn-]+//g;

					for my $i ( split(//, $rep_chars ) ){
						$line =~ s/$i/N/g;
					} 
				
					# feedback
					print "Warning: replacing unexpected characters ($rep_chars) with Ns for genome $genome\n";
				
					# store sequence
					push(@c_seq, $line);
				
				}
			
			}
		}
	
	}close INPUT;
	
	# Store final sequence.
	my $store_seq = join ("", @c_seq);
	$seq_store{$contig_id} = $store_seq;
	
	# Get co-ordinates and print sequences.
	open INPUT, "$gff_dir/$genome.gff" or die $!;
	while(<INPUT>){
	
		my $line = $_;
		chomp $line;
		my @line_array = split(/\t/, $line);
	
		# Variables.
		my $contig="";
		my $sta="";
		my $end="";
		my $strand="";
		my $id="";

		if( ($line !~ /^##/) && ($line !~ /^#!/) ){
		
			unless( $line_array[2] eq "gene"){
						
				# Set variables 
				$contig = $line_array[0];
				$sta = $line_array[3];
				$end = $line_array[4];
				my $type = $line_array[2];
		
				# Direction
				if($line_array[6] eq "+"){
					$strand="Forward";
				}elsif($line_array[6] eq "-"){
					$strand="Reverse";
				}
		
				# feature length
				my $len = (($end - $sta) + 1);
			
				# Clear variables
				$id = "";
			
				# Split info line.
				my @info = split (/;/, $line_array[8]);
				foreach( @info ){
				
					# Prokka
					if ($_ =~ /^locus_tag=(.+)/){
						$id = $1;
					}
					# RAST
					if ($_ =~ /^ID=(.+)/){
						$id = $1;
					}
					
				}
				
				# check if id is in an unstored gene family
				if ( ($id ne "") && ($loci_group { $id }) ){
				
					my @group_ids = keys %{$loci_group { $id }};
					
					for my $group_id (@group_ids){ 
					
						# store if representative has not already been stored or the stored sequence length is shorter
						if ( !$stored_groups{$group_id} || ( $stored_groups{$group_id} < $len ) ){
				
							# Get sequence from contig store.
							my $seq = substr( $seq_store{$contig}, $sta-1, $len );
					
							# Reverse complement if necessary.
							if ($strand eq "Reverse"){
								$seq = reverse($seq);
								$seq =~ tr/ATCG/TAGC/;
							}
						
							my $group_no = $group_list{$group_id};
							my $group_product = $products{ $group_id };
							my $group_gene = $gene_names{ $group_id };
						
							my $header = sprintf( ">%s;representative_genome=%s;locus_tag=%s;gene_name=%s;gene_product=%s;number_genomes=%s", $group_id, $genome, $id, $group_gene, $group_product, $group_no);
						
							# store info line and sequence
							$stored_seq {$group_id}{"s"} = $seq; 
							$stored_seq {$group_id}{"h"} = $header; 
						
							# Stored sequence added to sanity checking variable.
							$stored_groups {$group_id} = $len;
					
						}
					
					}
				
				}
		
			}
			
		}elsif($line =~ /^##FASTA/){
			last;
		}
		
	}close INPUT;
	
	# finish if all families found
	last if( scalar(keys(%stored_groups)) == scalar(keys(%group_list)) );
	 
}

# Check all sequences have been extracted.
for my $l_check ( sort keys %group_list ){
	print " - WARNING: No sequence found for gene_family/allele $l_check.\n" unless $stored_groups {$l_check};
}


# open output files
open OUTPUTN, ">$output_dir.ffn" or die " - ERROR: could not open output file\n";
open OUTPUTA, ">$output_dir.faa" or die " - ERROR: could not open output file\n";

# print sequence to file 
print " - printing sequences to files\n";
for my $k ( sort keys %stored_seq ) {

	my $seq = $stored_seq{$k}{"s"};
	my $header = $stored_seq{$k}{"h"};
	
	# convert to aa
	my $seq_obj = Bio::Seq->new( -seq => $seq, -alphabet => 'dna' );
	my $seq_aa = $seq_obj->translate()->seq();
	
	# print to files
	print OUTPUTN "$header\n$seq\n";
	print OUTPUTA "$header\n$seq_aa\n";
}

print " - complete\n";
	

exit
