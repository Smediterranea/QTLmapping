#!/usr/bin/perl
use strict;
use warnings;


#this code reads in a list of scaffolds (one scaffold per line) to convert it to a string separated by a space

#my $list = "f2snps_nomiss.gt.txt"; 
my $list = "f2.2parents.aabb.miss6.txt"; 


open (PRIMER, '<', $list)or die "can't read the genome file.\n";

my $i;

while (my $line = <PRIMER>){
	chomp $line;
	if ($line =~ /CHROM/){
		my @names = split "\t", $line;
		my @tmp;my $length = @names;
		for ($i=4; $i <$length; $i++){
			$names[$i] =~ /(^.+).GT/;
			$names[$i] = $1;
			push @tmp, $names[$i];
								}
		my $id = $names[0]."_".$names[1]."_".$names[2]."_".$names[3];
		my $space = "X";
		my $gt = join "\t", @tmp;
		print $id, "\t", $space, "\t", $gt, "\n";
							}
	else{
		my @all = split "\t", $line;
		my @b; my $length = @all;
		my $ref = $all[2]."/".$all[2];
		my $het = $all[2]."/".$all[3];
		my $alt = $all[3]."/".$all[3];
		for ($i=4; $i <$length; $i++){
			#if ($all[$i] eq $ref){$all[$i]=0;}
			#elsif ($all[$i] eq $alt){$all[$i]=1;}
			#elsif ($all[$i] eq $het){$all[$i]=2;}
			if ($all[$i] eq $ref){$all[$i]="A";}
			elsif ($all[$i] eq $alt){$all[$i]="B";}
			elsif ($all[$i] eq $het){$all[$i]="X";}
			else {$all[$i]="U";}
			push @b, $all[$i];
			}
		my $id = $all[0]."_".$all[1]."_".$all[2]."_".$all[3];
		$id =~ /dd_Smes_g4_(\d+)/;
		my $gt = join "\t", @b;
		print $id, "\t", $1, "\t", $gt, "\n";
		}
} 
