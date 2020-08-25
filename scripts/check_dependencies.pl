#!/usr/bin/env perl 

use strict;
use warnings;

# Check dependencies for PIRATE.
# Note: currently no version check.

# parallel
my $parallel = 0;
if( `command -v parallel;` ) {
    $parallel = 1;
} else {
    warn(" - ERROR: GNU parallel not found in system path.\n");
}

# cd-hit
my $cd_hit = 0;
for my $exe (qw(cdhit cd-hit) ) {
    if( `command -v $exe;` ) {
	$cd_hit = 1;
	last;
    }
}
warn(" - ERROR: cd-hit (or alternate invocation cdhit) not found in system path.\n") if $cd_hit == 0;

# blast+
my $blast = 0 ;
if ( `command -v blastn;` ) {
    $blast = 1;
} else {
    warn(" - ERROR: blast+ not found in system path.\n");
}

# mcl 
my $mcl = 0;
if( `command -v mcl;` ) {
    $mcl = 1;
} else {
    warn(" - ERROR: mcl not found in system path.\n");
}

# die if dependencies are not available.
unless ( $cd_hit && $blast && $mcl && $parallel ) { exit(2) }

# check optional dependencies

# diamond
unless ( `command -v diamond` ) {
    warn("\n - WARNING: cannot find diamond binary, cannot use --diamond command.\n");
}

# R

unless( `command -v R;`) {
    warn(" - WARNING: R not found in system path, cannot use -r command.\n");
}

# fasttree
my $ft = 0;
for my $exe (qw(FastTreeMP FastTree fasttree) ) { 
    if( `command -v $exe;`) {
	$ft = 1;
	last;
    }
}

warn(" - WARNING: fasttree not found in system path, a binary presence-absence tree will not be created.\n") if $ft == 0;
exit
