#!/usr/bin/perl

use strict;
use warnings;


#this code reads in a list of scaffolds (one scaffold per line) to convert it to a string separated by a space

my $list = "parents.gt.txt"; 

open (PRIMER, '<', $list)or die "can't read the genome file.\n";

my $i;

#output AA/AB genotypes
while (my $line = <PRIMER>){
	chomp $line;
	if ($line =~ /CHROM/){print $line, "\n";}
	else{
		my @all = split "\t", $line;
		my $ref = $all[2]."/".$all[2];
		my $het = $all[2]."/".$all[3];
		my $alt = $all[3]."/".$all[3];
		if ($all[4] eq $ref and $all[5] eq $het) {print $line,"\n";}
				elsif ($all[4] eq $alt and $all[5] eq $het) {print $line,"\n";}
		}
} 






__END__
#output AA/BB genotypes
while (my $line = <PRIMER>){
	chomp $line;
	if ($line =~ /CHROM/){print $line, "\n";}
	else{
		my @all = split "\t", $line;
		my $ref = $all[2]."/".$all[2];
		my $het = $all[2]."/".$all[3];
		my $alt = $all[3]."/".$all[3];
		if ($all[4] eq $het or $all[4] eq "./.") {next;}
				elsif ($all[5] eq $het or $all[5] eq "./.") {next;}
				elsif ($all[5] eq $all[4]) {next;}
				else {print $line,"\n";}

		}
} 
