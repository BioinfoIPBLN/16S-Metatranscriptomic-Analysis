#!/usr/bin/perl

use strict;
use Getopt::Long;
use Sort::Key::Natural qw( natsort );
my $dir;
my $help;
my $label;

GetOptions(
        "help" => \$help,
        "dir|d=s" => \$dir,
		"label|l=s" => \$label
);
if($dir and $label){
        #Declaration of the output path 
        my $fileresults= $label. "_counts.tab"; 
        my $fileresultsOTU= $label. "_OTU.tsv"; 
        #Variable declaration
        my $hashref;
        my $data;
        my $files;
        my $taxa_info;
        print STDERR date(). " :: Reading $dir\n";
        # Opening and reading the directory supplied by the user
		my @files=split(" ",$dir);
        # #Reading each file of the directory
        foreach my $file( sort @files){
             if ($file =~ "bracken"){
                #Opening each file of the directory
                open (FILE, $file) || die "Can't open '$file': $!";
                #Saving the files name in a hash
                #print "$file\n";
			 	my ($sample,undef)=split("_bracken",$file);
				$sample=~s/.*\/(.*)/$1/g;
				#print STDERR $sample ."\n";
                #print $sample ."\n";
			 	$files->{$sample}++;
                #Reading the file. Each line contains the number of the miRNA, a tab and the number of reads
                while(<FILE>){
                    #Deleting the last character of the line
                    chomp;
                    if($_ !~ /^name/){
                        #Using the split with the tab we can keep the name of the miRNA in $miRNA variable
                        #and the number of reads in $reads variable. 
                        my @data=split(/\t/);
                        my $miRNAs=$data[0];
                        $miRNAs=~s/ /_/g;
                        $miRNAs=~s/\'//g;
                        $miRNAs=~s/\#//g;

                        my $reads=$data[5];
                        my $taxa=$data[2];

                        $hashref->{$miRNAs}++;
                        $data->{$miRNAs}->{$sample}="$reads";   
                        $taxa_info->{$miRNAs}=$taxa;
			 			
			 			#print "$miRNAs\t$sample\t$reads\n";
                    }
                }
                #Closing the file
                close FILE;
             }
        }
        close DIR;
        print STDERR date(). " :: Creating $fileresults\n";
        print STDERR date(). " :: Creating $fileresultsOTU\n";
        #Opening the output file to write the data
        open(RESULTS,"> ".$fileresults) || die $!;
        open(OTU,"> ".$fileresultsOTU) || die $!;
        #Printing on the first row of the results file the name of the input files sorted and 
        #separated by a tab 
        print RESULTS "OTU\t". join("\t",(natsort keys (%$files)))."\n";
        print OTU "OTU\tTaxonomy\n";
        #The next lines will be printed going over the hash. Reading the first dimension of the hash
        my $cont=1;
        foreach my $miRNAs(sort keys %$data){
                #Printing the name of each miRNA
                print RESULTS $cont."\t";
                print OTU "$cont\t$miRNAs\t".$taxa_info->{$miRNAs} ."\n";
                my @results;
                #An array is used to save the number of reads of each file for a defined miRNA.  
                foreach my $file(natsort keys %{$files}){
                        my $valor=0;
                        if(exists $data->{$miRNAs}->{$file}){
                                $valor=$data->{$miRNAs}->{$file};
                        }
                        push(@results,$valor);
                }
                #Printing the number of the reads of each miRNA separated by a tab 
                print RESULTS join("\t",@results)."\n";
                $cont++;
        }
        #Closing the file
        close RESULTS;
        close OUT;
        print STDERR date(). " :: Finished, check $fileresults\n";

}
elsif($help){
        help();
}
else{
        help();
}

sub help{
        print STDERR "\n\tperl $0 -dir=ESCA-ReadCoun/ -label ESCA\n\n";
}

sub date{
        my $date=("date \"+%D %H:%M:%S\"");
        $date=`$date`;
        chomp $date;
        return("[".$date."]");
}
