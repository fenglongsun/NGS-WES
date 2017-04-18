#! /usr/bin/perl -w
@hs=(6,7,8,10,15,16,17,19,20,22,23,24,25,27,28,29,30,31,32,33,35,36,37,38,39,40,41,42);
@alig= (1,2,3,4,7,9,13,14,15,18,20,21,22);
open IN,$ARGV[0]||die;
while(<IN>){
	chomp;
	if(/^BAIT/){
		@tt=split(/\s+/);
	}elsif(/^Bait/){
		@ss=split(/\s+/);	
	}elsif(/^CATE/){
		@al1=split(/\s+/);
	}elsif(/^PAIR/){
		@al2=split(/\s+/);
	}elsif(/TOTAL_READS/){
		@qc1=split(/\s+/);
	}elsif(/\d+/){
		@qc2=split(/\s+/);
	}
}
print "Sample_ID\t$ARGV[2]\n";
if($ARGV[1] eq 'hs'){
foreach(@hs){
	print "$tt[$_ - 1]\t$ss[$_ - 1]\n";
}}
if($ARGV[1] eq 'alig'){
foreach(@alig){
	print "$al1[$_ - 1]\t$al2[$_ - 1]\n";
}}
if($ARGV[1] eq 'qc'){
foreach $i(0..10){
	
        print "$qc1[$i]\t$qc2[$i]\n";
}
}
