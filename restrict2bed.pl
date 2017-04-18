#! /usr/bin/perl -w
#perl restrict2bed.pl [.bed] [.vcf] > [restricted.vcf]
open IN,$ARGV[0]||die;
while(<IN>){
	chomp;
	@tt=split/\s+/;
	push @{$info{$tt[0]}},[$tt[1],$tt[2]];
}
my %out;
open IN1,$ARGV[1]||die;
while(<IN1>){
	chomp;
	if(/^#/){print "$_\n";next};
	$line=$_;
	if($line!~/^chr/){
        	$line="chr".$line;
        }
        @tt=split(/\s+/,$line);
	$key=$tt[0];
	$kk=$tt[0]."\t".$tt[1];	
	if(defined $info{$key}){
		@ss=@{$info{$key}};
	}else{
		next;
	}
	foreach my $i(0..$#ss){
		if($tt[1]>=$ss[$i][0] && $tt[1]<=$ss[$i][1] && (! defined $out{$kk})){
			push @{$xx{$key}},[$tt[1], $line];
			$out{$kk}++;
		}
	}
}
foreach my $chr( sort keys %xx){
	@mm=@{$xx{$chr}};
	@mm=sort {$a->[0] <=> $b->[1]} @mm;
	foreach my $i(0..$#mm){
		print "$mm[$i][1]\n";
	}
}
