#!usr/bin/perl -w
my $input=$ARGV[0];
open (IN,$input) or die "$!";
my %hash;my %ch;
while (<IN>){
	if ($_=~/sample/){
		print "$_";	
		next;
	};
	my $line=$_;
	my @tmp=split /\t+/,$_;
	my $sample=$tmp[0];
	my $chr=$tmp[1];
	my $start=$tmp[2];
	my $end=$tmp[3];
	next if ($end-$start<600000); 
	#my $num=$tmp[4];
	my $mean=$tmp[4];
	$hash{$sample}{$chr}{$start}=$line;
	$ch{$sample}{$chr}++;
}
#open (OUT,">sample_chr.list") or die "$!";
#foreach my $sample (sort keys %ch){
#	my @tmp;
#	foreach my $chr (sort keys %{$ch{$sample}}){
#		push @tmp,$chr;
#	} 
#	my $num=@tmp;
#	print "$sample\t$num\n";
#}
my @chrs=(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22);
foreach my $sample(sort keys %hash){
	for my $i (0,1..$#chrs){
		if ($ch{$sample}{$chrs[$i]}){
			foreach my $start (sort{$a<=>$b} keys %{$hash{$sample}{$chrs[$i]}}){
				my $line=$hash{$sample}{$chrs[$i]}{$start};
				my ($sample,$chr,$start,$end)=split /\t+/,$line;
				if (($end - $start)>0){
					print "$line";
				}
				else{
					print "NA\n";
				}
			}
		}
		else {
			print "NA\n";
		}
		
	}
}
#print join ("\t",@chrs)."\n";
