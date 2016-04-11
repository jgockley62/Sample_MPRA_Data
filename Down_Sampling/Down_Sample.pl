#usr/bin/perl 

use warnings;
use strict;
use POSIX;

#This script will down sample the data in the MPRA sequencing file proportionatly

# an array will be created where column six is the number of occurences of column 1
#
#USAGE: perl Down_Sample.pl Rep2_mRNA.txt
#Input File: Name       Hit Sequence    Tag-Sequence    Observed_Counts Proportion_of_Total_Counts      Whole_Number_Proportion
#Input File: synCRE_Promega_0   GCACCAGACAGTGACGTCAGCTGCCAGATCCCATGGCCGTCATACTGTGACGTCTTTCAGACACCCCATTGACGTCAATGGGAGAAC ATCTGCGGCC      25      2.1629871745762e-07     2163

my $hand1=$ARGV[0];

#Total counts in full scale experiment
my $TotalCounts=0;

#Same as input but addl columns will have downsampled counts
my %Final=();

#Taglist to print out final Table in order
my @FragList=();

#Array to permute through
my @Permer=();

open FILE1, "$hand1" or die "Could not open $hand1: $!\n";
while(my $tmp = <FILE1>){
    $tmp =~ s/[\n\r]//g;
    my @line = split("\t", $tmp);
    
    #Add to total Observed Counts
    $TotalCounts = $TotalCounts+$line[3];
    
    #Create entry for final table
    $Final{$line[0]} = $tmp."\t0\t0\t0\t0\t0";
        
    #Add FragList Entry
    push(@FragList,$line[0]);
    
    #Seed permution template
    for (my $i = 0; $i<$line[5]; $i++){
        push(@Permer,$line[0]);
    }
}
close FILE1;
print "File Loaded\n";
#Number of options in probability array
my $SIZE = scalar @Permer;
print $SIZE."\n";
print localtime."\n";

#DownSample to 5%
for (my $i = 0; $i<ceil(0.05*$TotalCounts); $i++){
        
        my $Tag = $Permer[ int(rand($SIZE))-1 ];
        my @temp = split("\t", $Final{$Tag});
        #Add a read count
        $temp[10] = $temp[10]+1;    
        #my $FOOs = join "\t", @temp;
        $Final{$Tag} = join "\t", @temp;

    }
print "5% Perms Done\n"; 
print localtime."\n";
#DownSample to 10%
for (my $i = 0; $i<ceil(0.1*$TotalCounts); $i++){
        
        my $Tag = $Permer[ int(rand($SIZE))-1 ];
        my @temp = split("\t", $Final{$Tag});        
        #Add a read count
        $temp[9] = $temp[9]+1; 
        $Final{$Tag} = join "\t", @temp;
    }
print "10% Perms Done\n";
print localtime."\n";
#DownSample to 25%
for (my $i = 0; $i<ceil(0.25*$TotalCounts); $i++){
        
        my $Tag = $Permer[ int(rand($SIZE))-1 ];
        my @temp = split("\t", $Final{$Tag});        
        #Add a read count
        $temp[8] = $temp[8]+1;
        $Final{$Tag} = join "\t", @temp;
    }
print "25% Perms Done\n";
print localtime."\n";
#DownSample to 50%
for (my $i = 0; $i<ceil(0.5*$TotalCounts); $i++){
        
        my $Tag = $Permer[ int(rand($SIZE))-1 ];
        my @temp = split("\t", $Final{$Tag});        
        #Add a read count
        $temp[7] = $temp[7]+1; 
        $Final{$Tag} = join "\t", @temp;   
    }
print "50% Perms Done\n";  
print localtime."\n";
#DownSample to 75%
for (my $i = 0; $i<ceil(0.75*$TotalCounts); $i++){
        
        my $Tag = $Permer[ int(rand($SIZE))-1 ];
        my @temp = split("\t", $Final{$Tag});
        #Add a read count
        $temp[6] = $temp[6]+1;
        $Final{$Tag} = join "\t", @temp;
    }
print "75% Perms Done\n";
print localtime."\n"; 

#####Print Out Permuted Data
##Create Output File
my @Name = split(/\./, $hand1);
my $OUT = $Name[0]."_DownSampled.txt";
open(my $out, '>', $OUT) or die "Could not open file $OUT $!\n";
#Print out final Matrix
print $out "Name\tHit_Sequence\tTag_Sequence\tObserved_Counts\tProportion_of_Total_Counts\tWhole_Number_Proportion\t75%\t50%\t25%\t10%\t5%\n";
foreach my $Fragment (@FragList){

    print $out $Final{$Fragment}."\n";

}

