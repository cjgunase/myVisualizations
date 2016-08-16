#use strict;
use warnings;
use List::Util qw( min max );
use POSIX;
my(%chrLenHash);
open(IN,"chr_len.txt") or die $!;
while(<IN>){
	my @hashVal = split('\t', $_);
	$chrLenHash{$hashVal[0]}=$hashVal[1];
}
close(IN);

@fileArray=("birch_chr.gff","birch_mite.gff","copia_chr.gff","DNA_type_chr.gff","gypsy_chr.gff","line_chr.gff","sine_chr.gff","SSR_chr.gff");

foreach my $key (@fileArray){

my $LGNUM =14;
#qx(rm "birch_heat.txt");

		for(my $i=1;$i<=$LGNUM;$i++){
		my @array=();
		open(IN,$key) or die " file cannot open$!";#change file here
		open(OUT,">tempLG$i.txt") or die $!;
		while(<IN>){
			s/\r\n?/\n/;  #remove return
			chomp();
			 my @values = split('\t', $_);
			 if($values[0] eq "LG$i"){
				print OUT "$_\n";
			 }
		 
		}
		
	close(OUT);
	close(IN);
	qx(sort -nk 4,4 tempLG$i.txt>LG$i.txt);#sort the whole this by ascending values
	qx(rm tempLG$i.txt);	
	
	$stepsize=1000000;
	my $n=floor($chrLenHash{"LG$i"}/$stepsize);
	my $range=0;
	my $count=0;
	my $point=0;
	my $j=0;
	
	open(OUT,">>$key"."_heat.txt") or die $!;
	while($j < $n){

	$range+=$stepsize;

	open(IN,"LG$i.txt") or die " file cannot open$!";
	$file=0;
	while(<IN>){
	$file++;
	s/\r\n?/\n/;  #remove return
	chomp();	
	unless($file <= $point){	
		my @values = split('\t', $_);
						
			if($values[2] eq "Transposon" or  $values[2] eq "TEprotein" or $values[2] eq "CDS" or $values[1] eq "mite"){
				if($values[4] < $range){
				$count++;
				$point++;
				}else{
				last;
				}
			 
			}
	}	
		 
		 
	}

	$j++;
	#print OUT "LG$i\t$count\n";
	push @array,$count;
	$count=0;
	close(IN);

	}
	$min = min @array;
	$max = max @array;
	
	
	$start=0;
	$end=49999;
	foreach(@array){
	
	$heat=($_-$min)/($max-$min);
	print OUT "LG$i\t$start\t$end\t$heat\n";
	$start=$end+1;
	$end = $end+$stepsize;
	
	}
	qx(rm LG$i.txt);

	}
	
	

	
	}
	