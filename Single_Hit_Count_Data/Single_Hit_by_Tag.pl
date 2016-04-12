#usr/bin/perl 

use warnings;
use strict;
use POSIX;

#This script will reformat data to one Row Per TAG


#USAGE: perl Single_Hit_bt_Tag.pl SingleHit_CountData_RAW.txt

my $hand1=$ARGV[0];

#Replicate 1
my $OUT = "Tags_Replicate_1.txt";
open(my $Rep1, '>', $OUT) or die "Could not open file $OUT $!\n";
#Header
print $Rep1 "Name\tTag\tmRNA_Counts\tpDNA_Counts\n";

#Replicate 1
my $OUT2 = "Tags_Replicate_2.txt";
open(my $Rep2, '>', $OUT2) or die "Could not open file $OUT $!\n";
#Header
print $Rep2 "Name\tTag\tmRNA_Counts\tpDNA_Counts\n";

open FILE1, "$hand1" or die "Could not open $hand1: $!\n";
while(my $tmp = <FILE1>){
    $tmp =~ s/[\n\r]//g;
    my @line = split("\t", $tmp);
    if($line[0] eq 'Variant_ID/Tag_ID'){
    }else{
        #Seed permution template
        for (my $i = 1; $i<14; $i++){
            print $Rep1 $line[0]."\t".$i."\t".$line[$i]."\t".$line[$i+13]."\n";
            print $Rep2 $line[0]."\t".$i."\t".$line[$i+26]."\t".$line[$i+39]."\n";
        }
    }
}
