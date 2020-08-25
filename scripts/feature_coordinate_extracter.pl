#!/usr/bin/env perl

# Extract feature co-ordinates from GFF file.

# Dependencies
use strict; 
use warnings;
use Getopt::Long qw(GetOptions :config no_ignore_case);
use Pod::Usage;
use Cwd 'abs_path';
use File::Basename;

# Version

=head1  SYNOPSIS


=head1 Descriptions
	
...

=cut

# command line options
my $man = 0;
my $help = 0;
my $input_file = '';
my $output_file = '';
my $output_intergenic = '';

my $feature_list = "CDS";

GetOptions(
	'help|?' 	=> \$help,
	'man' 		=> \$man,
	'input=s' 	=> \$input_file,
	'output=s'	=> \$output_file,
	'features=s'    => \$feature_list,
	'intergenic=s' => \$output_intergenic,
) or pod2usage(2);
pod2usage(1) if $help;
pod2usage(-verbose => 2) if $man;

# Make regex for feature test.
my @features = split(/,/, $feature_list);
my $regex = sprintf("\(%s\)", join("|", @features));

# check file exists
die "ERROR: input file not found.\n" unless -f $input_file;

# isolate named after input GFF
my $isolate = basename($input_file,(".modified.gff", ".gff"));

# Outputs
open OUTPUT_G, ">$output_file" or die "ERROR: Could not open $output_file.\n";
print OUTPUT_G join("\t",qw(Name Gene Start End Length Type Strand Contig Product)),"\n";

# Optionally produce intergenic co-ordinates.
unless( $output_intergenic eq '' ){
	open OUTPUT_I, ">$output_intergenic" or die "ERROR: Could not open $output_intergenic.\n";
	print OUTPUT_I join("\t",qw(Name Gene_name Start End Length Type Contig Product)),"\n";
}

# Parse input
my @gene_array = ();
my %contig_hash_end = ();
open INPUT, "$input_file" or die "ERROR: Could not open input file - $input_file\n";
while(my $line = <INPUT>){
	chomp $line;
	
	my @line_array = split(/\t/, $line);
	
	# variables
	my $sta="";
	my $end="";
	my $strand="";
	my $id="";
	my $gene="";
	my $product = "";
	my $type = "";

	my $contig="";
	my $seq_end = "";
		
	if( ($line !~ /^##/) && ($line !~ /^#!/) ) {
		if( $line_array[2] eq "gene"){
			# ignore genes
		} elsif ($line_array[2] =~ /^$regex/){
			my %group = map { split(/=/, $_, 2) } split(/;/,$line_array[8]);
			($contig,undef,$type,$sta,$end) = @line_array;
		
			$strand = ($line_array[6] eq "+") ? "Forward" : "Reverse";
		
			my $len = (($end - $sta) + 1);
			if( $group{ID} =~ /${isolate}/ ) {
				$id = $group{ID}
			}

			$gene = ( exists $group{gene} ) ? $group{gene} : "";
			if ( exists $group{product} ) {
				$product = $group{product};
			} elsif( exists $group{Name} ) {
				$product = $group{Name};
			} else {
				$product = $group{Parent};
			}
			# RAST gffs
			$gene = $group{Parent} if ( ! length($gene) && exists $group{Parent} );
			# store feature co-ordinates.
			my @tmp_array=("$id", "$gene", "$sta", "$end", "$len", "$type", "$strand", "$contig");
			push @gene_array, [@tmp_array];
			
			# print output
			print OUTPUT_G join("\t",$id,$gene,$sta,$end,$len,$type,$strand,$contig,$product),"\n";
		}
	} elsif($line =~ /^##sequence-region\s+(\S+)\s+(\d+)\s+(\d+)/){
		($contig,$seq_end) = ($1,$3);
		$contig_hash_end{$contig} = $seq_end;
	} elsif($line =~ /^##FASTA/){
		last;
	}
}

# if only genic co-ordinates are needed stop here.
unless ( $output_intergenic eq '' ){ exit }

# Otherwise indentify intergenic co-ordinates. 
#$gene_count=scalar(@gene_array);

#$int_count=0;
#for($i=0; $i<$gene_count; $i++){
#	# First gene.
#	if($i == 0){
#		$contig=$gene_array[$i][7];
#		$int_sta=1;
#		$int_end=($gene_array[$i][2] - 1);
#		$int_len=(($int_end - $int_sta) + 1);
#		
#		$gene_pre="NA";
#		$gene_pos=$gene_array[$i][0];
#		$int_type="NA";
#		$int_name="${gene_pre}_+_+_${gene_pos}_+_+_$int_type";
#		
#		$int_count++;
#		$int_id="${isolate}_intergenic_$int_count";
#		
#		if($int_len > 0){
#			print OUTPUT_I "$int_id\t$int_name\t$int_sta\t$int_end\t$int_len\t$int_type\t$contig\tNA\n";
#		}
#	}
#	# Last gene.
#	elsif($i == ($gene_count - 1)){
#		$contig=$gene_array[$i][7];
#		$int_sta=($gene_array[$i][3] + 1);
#		$int_end=$contig_hash_end{$contig};
#		$int_len=(($int_end - $int_sta) + 1);
#		
#		$gene_pre=$gene_array[$i][0];
#		$gene_pos="NA";
#		$int_type="NA";
#		$int_name="${gene_pre}_+_+_${gene_pos}_+_+_$int_type";
#		
#		$int_count++;
#		$int_id="${isolate}_intergenic_$int_count";
#		
#		if($int_len > 0){
#			print OUTPUT_I "$int_id\t$int_name\t$int_sta\t$int_end\t$int_len\t$int_type\t$contig\tNA\n";
#		}
#	}
#	# Gene on different contig.
#	elsif($gene_array[($i-1)][7] ne $gene_array[$i][7]){
#		# Edge of previous contig.
#		$contig=$gene_array[($i-1)][7];
#		$int_sta=($gene_array[($i-1)][3] + 1);
#		$int_end=$contig_hash_end{$contig};
#		$int_len=(($int_end - $int_sta) + 1);
#		
#		$gene_pre=$gene_array[($i-1)][0];
#		$gene_pos="NA";
#		$int_type="NA";
#		$int_name="${gene_pre}_+_+_${gene_pos}_+_+_$int_type";
#		
#		$int_count++;
#		$int_id="${isolate}_intergenic_$int_count";
#		
#		if($int_len > 0){
#			print OUTPUT_I "$int_id\t$int_name\t$int_sta\t$int_end\t$int_len\t$int_type\t$contig\tNA\n";
#		}
#		
#		# Edge of next contig.
#		$contig=$gene_array[$i][7];
#		$int_sta=1;
#		$int_end=($gene_array[$i][2] - 1);
#		$int_len=(($int_end - $int_sta) + 1);
#		
#		$gene_pre="NA";
#		$gene_pos=$gene_array[$i][0];
#		$int_type="NA";
#		$int_name="${gene_pre}_+_+_${gene_pos}_+_+_$int_type";
#		
#		$int_count++;
#		$int_id="${isolate}_intergenic_$int_count";
#		
#		if($int_len > 0){
#			print OUTPUT_I "$int_id\t$int_name\t$int_sta\t$int_end\t$int_len\t$int_type\t$contig\tNA\n";
#		}
#	}
#	# Gene on same contig.
#	elsif($gene_array[($i-1)][7] eq $gene_array[$i][7]){
#		$contig=$gene_array[$i][7];
#		$int_sta=($gene_array[($i-1)][3] + 1);
#		$int_end=($gene_array[$i][2] - 1);
#		$int_len=(($int_end - $int_sta) + 1);
#		
#		$gene_pre=$gene_array[($i-1)][0];
#		$gene_pos=$gene_array[$i][0];
#	
#		if($gene_array[($i-1)][6] eq "Forward" && $gene_array[$i][6] eq "Forward"){
#			
#			$int_type="CO_F";
#			$int_name="${gene_pre}_+_+_${gene_pos}_+_+_$int_type";
#		}elsif($gene_array[($i-1)][6] eq "Forward" && $gene_array[$i][6] eq "Reverse"){
#			
#			$int_type="DT";
#			$int_name="${gene_pre}_+_+_${gene_pos}_+_+_$int_type";
#		}elsif($gene_array[($i-1)][6] eq "Reverse" && $gene_array[$i][6] eq "Forward"){
#			
#			$int_type="DP";
#			$int_name="${gene_pre}_+_+_${gene_pos}_+_+_$int_type";
#		}elsif($gene_array[($i-1)][6] eq "Reverse" && $gene_array[$i][6] eq "Reverse"){
#			
#			$int_type="CO_R";
#			$int_name="${gene_pre}_+_+_${gene_pos}_+_+_$int_type";
#		}
#		
#		$int_count++;
#		$int_id="${isolate}_intergenic_$int_count";
#		
#		if($int_len > 0){
#			print OUTPUT_I "$int_id\t$int_name\t$int_sta\t$int_end\t$int_len\t$int_type\t$contig\tNA\n";
#		}
#	}
#}
#
exit

