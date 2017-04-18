#modified by zhangyu, 20160217
#!/usr/bin/perl -w
use strict;
use lib '/usr/local/share/perl5/';
use lib '/home/zhanghk/software/perl/lib/share/perl5';
use lib '/opt/perl/lib/site_perl/5.14.2';
use Statistics::R;
use Getopt::Long;
use File::Basename;

my ($annovar, $stat, $type, $help);

GetOptions(
	'a=s'  => \$annovar,
	's=s' => \$stat,
	't=s' => \$type,
	'help' => \$help
);

my $usage=<<INFO;
Usage:
    perl $0 [options]
Options:
           -a   <file>    input annovar annotation result, SampleName_.*multianno.txt
           -s   <file>    output file name of statistic information
           -t   <string>  the type of variation, include(SNV/SNP/INDEL)
Example:
	perl annovar_statistic.pl -a cnv.txt -s cnv.stat -t CNV
Authro:
	zhanghaikuan\@berrygenomics.cn
INFO

if($help || !$annovar || !$stat){
	die $usage;
}
my $dir = dirname $annovar;
my $sample_name = $stat;
my (%func_hash,%draw_snp) = ();
my ($total,$ti,$tv) = (0,0,0);
my (@indel,@insertion,@deletion);

open IN,"<$annovar" or die "input annovar result open error! $!\n";
my $head = <IN>;
$head =~ /avsnp(\d+).*/;
my $dbsnp_version = "dbsnp$1";
while(my $line = <IN>){
	chomp $line;
	my @fields = split /\t+/,$line;
	my ($chr, $start, $end, $ref, $alt, $func, $gene, $exonic, $g1000, $dbsnp) = ($fields[0], $fields[1], $fields[2], $fields[3], $fields[4], $fields[5], $fields[6], $fields[8], $fields[18], $fields[11]);

	if(($type =~ /SNV/i || $type =~ /SNP/i)&&($ref ne "-" && $alt ne "-")){	
		$func_hash{$func}++;
		if($g1000 ne "." and $dbsnp ne "."){$func_hash{"1000genome and $dbsnp_version"}++;}
		if($g1000 ne "." and $dbsnp eq "."){$func_hash{"1000genome only"}++;}
		if($g1000 eq "." and $dbsnp ne "."){$func_hash{"$dbsnp_version only"}++;}
		if($dbsnp eq "." and $g1000 eq "."){$func_hash{"Novel"}++;}
		if($exonic ne "."){
			$func_hash{$exonic}++;
		}
		$total++;
		if($type =~ /SNP/i || $type =~ /INDEL/i){if($line =~ /0\/0\:?/ or $line =~ /1\/1:?/){$func_hash{"Homo"}++;}else{$func_hash{"Het"}++;}}
		if($type =~ /SNV/i){if($line =~ /\thom\t/){$func_hash{"Homo"}++;}else{$func_hash{"Het"}++;}}


		if((($ref eq "A") and ($alt eq "G")) or (($ref eq "G") and ($alt eq "A")) or (($ref eq "C") and ($alt eq "T")) or (($ref eq "T") and ($alt eq "C"))){
			$ti++;
		}else{
			$tv++;
		}
		if((($ref eq "T") and ($alt eq "C")) or (($ref eq "A") and ($alt eq "G"))){
		$draw_snp{"T:A->C:G"}++;
		}
		if((($ref eq "T") and ($alt eq "A")) or (($ref eq "A") and ($alt eq "T"))){
			$draw_snp{"T:A->A:T"}++;
		}
		if((($ref eq "T") and ($alt eq "G")) or (($ref eq "A") and ($alt eq "C"))){
			$draw_snp{"T:A->G:C"}++;
		}
		if((($ref eq "C") and ($alt eq "T")) or (($ref eq "G") and ($alt eq "A"))){
			$draw_snp{"C:G->T:A"}++;
		}
		if((($ref eq "C") and ($alt eq "G")) or (($ref eq "G") and ($alt eq "C"))){
			$draw_snp{"C:G->G:C"}++;
		}
		if((($ref eq "C") and ($alt eq "A")) or (($ref eq "G") and ($alt eq "T"))){
			$draw_snp{"C:G->A:T"}++;
		}
	}
	if($type =~ /INDEL/i&&($ref eq "-" || $alt eq "-")){
		$func_hash{$func}++;
		if($g1000 ne "." and $dbsnp ne "."){$func_hash{"1000genome and $dbsnp_version"}++;}
		if($g1000 ne "." and $dbsnp eq "."){$func_hash{"1000genome only"}++;}
		if($g1000 eq "." and $dbsnp ne "."){$func_hash{"$dbsnp_version only"}++;}
		if($dbsnp eq "." and $g1000 eq "."){$func_hash{"Novel"}++;}
		if($exonic ne "."){
			$func_hash{$exonic}++;
		}
		$total++;
		if($type =~ /SNP/i || $type =~ /INDEL/i){if($line =~ /0\/0\:?/ or $line =~ /1\/1:?/){$func_hash{"Homo"}++;}else{$func_hash{"Het"}++;}}
		if($type =~ /SNV/i){if($line =~ /\thom\t/){$func_hash{"Homo"}++;}else{$func_hash{"Het"}++;}}
		if($ref eq "-"){
			my $len = length($alt);
			$insertion[$len]++;
			$indel[$len]++;
		}elsif($alt eq '-'){
			my $len = length($ref);
			$deletion[$len]++;
			$indel[$len]++;
		}else{
			my $len1=length($ref);
			my $len2=length($alt);
			if($len1 < $len2){
				my $len=$len2-$len1;
				$insertion[$len]++;
				$indel[$len]++;
			}else{
				my $len=$len1-$len2;
				$deletion[$len]++;
				$indel[$len]++;
			}
		}
	}
}
close IN;

open OUT,">$dir/$stat\.$type\_statistic.xls" or die "STAT create error! $!\n";
print OUT "Total\t$total\n";
if($type =~ /SNV/i || $type =~ /SNP/i || $type =~ /INDEL/i){
	if(exists $func_hash{"1000genome and $dbsnp_version"}){
		my $rate = 100*$func_hash{"1000genome and $dbsnp_version"}/$total;
		printf OUT "%s\t%d(%4.2f%%)\n", "1000genome and $dbsnp_version", $func_hash{"1000genome and $dbsnp_version"}, $rate;
		delete $func_hash{"1000genome and $dbsnp_version"};
	}
	if(exists $func_hash{"1000genome only"}){
		my $rate = 100*$func_hash{"1000genome only"}/$total;
		printf OUT "%s\t%d(%4.2f%%)\n", "1000genome only", $func_hash{"1000genome only"}, $rate;
		delete $func_hash{"1000genome only"};
	}
	if(exists $func_hash{"$dbsnp_version only"}){
		my $rate = 100*$func_hash{"$dbsnp_version only"}/$total;
		printf OUT "%s\t%d(%4.2f%%)\n", "$dbsnp_version only", $func_hash{"$dbsnp_version only"}, $rate;
		delete $func_hash{"$dbsnp_version only"};
	}
	if(exists $func_hash{"Novel"}){
		my $rate = 100*$func_hash{"Novel"}/$total;
		printf OUT "%s\t%d(%4.2f%%)\n", "Novel", $func_hash{"Novel"}, $rate;
		delete $func_hash{"Novel"};
	}
	if(exists $func_hash{"Het"}){
		my $rate = 100*$func_hash{"Het"}/$total;
		printf OUT "%s\t%d(%4.2f%%)\n", "Het", $func_hash{"Het"}, $rate;
		delete $func_hash{"Het"};
	}
	if(exists $func_hash{"Homo"}){
		my $rate = 100*$func_hash{"Homo"}/$total;
		printf OUT "%s\t%d(%4.2f%%)\n", "Homo", $func_hash{"Homo"}, $rate;
		delete $func_hash{"Homo"};
	}
	if(exists $func_hash{"synonymous SNV"}){
		my $rate = 100*$func_hash{"synonymous SNV"}/$total;
		printf OUT "%s\t%d(%4.2f%%)\n", "synonymous SNV", $func_hash{"synonymous SNV"}, $rate;
		delete $func_hash{"synonymous SNV"};
	}
	if(exists $func_hash{"nonsynonymous SNV"}){
		my $rate = 100*$func_hash{"nonsynonymous SNV"}/$total;
		printf OUT "%s\t%d(%4.2f%%)\n", "nonsynonymous SNV", $func_hash{"nonsynonymous SNV"}, $rate;
		delete $func_hash{"nonsynonymous SNV"};
	}
	if(exists $func_hash{"exonic"}){
		my $rate = 100*$func_hash{"exonic"}/$total;
		printf OUT "%s\t%d(%4.2f%%)\n", "exonic", $func_hash{"exonic"}, $rate;
		delete $func_hash{"exonic"};
	}
	if(exists $func_hash{"intergenic"}){
		my $rate = 100*$func_hash{"intergenic"}/$total;
		printf OUT "%s\t%d(%4.2f%%)\n", "intergenic", $func_hash{"intergenic"}, $rate;
		delete $func_hash{"intergenic"};
	}
	if(exists $func_hash{"intronic"}){
		my $rate = 100*$func_hash{"intronic"}/$total;
		printf OUT "%s\t%d(%4.2f%%)\n", "intronic", $func_hash{"intronic"}, $rate;
		delete $func_hash{"intronic"};
	}
	
	foreach my $key (sort keys %func_hash){
		print OUT "$key\t$func_hash{$key}\n";
	}
	if($type !~ /INDEL/i){
		my $rate_ti_tv = $ti/$tv;
		printf OUT "Ti/Tv\t%4.2f\n",$rate_ti_tv;
	}
}
close OUT;

if($type =~ /SNP/i || $type =~ /SNV/i){
	#my $dir = dirname $stat;
	my $picture_snv = $dir."/".$sample_name."\.$type\_Spectrum.png";
	my $R = Statistics::R->new();
	$R->set( 'number', [$draw_snp{"T:A->C:G"},$draw_snp{"T:A->A:T"},$draw_snp{"T:A->G:C"}, $draw_snp{"C:G->T:A"},$draw_snp{"C:G->G:C"},$draw_snp{"C:G->A:T"}] );
	$R->set( 'name', ["T:A->C:G", "T:A->A:T", "T:A->G:C", "C:G->T:A", "C:G->G:C", "C:G->A:T"]);
	$R->run(qq`options(bitmapType='cairo')`); 
	$R->run(qq`png("$picture_snv", width=700,heigh=700)`);
	$R->run(qq`bar=barplot(number,ylim=c(0,1.35*max(number)),names.arg=name,cex=1.3,xlab="Mutation type",ylab="Mutation number",main="Mutation Spectrum",col=c(rgb(186,85,211,max=255),rgb(65,105,225,max=255),rgb(105,105,105,max=255),rgb(154,205,50,max=255),rgb(0,139,139,max=255),rgb(218,165,32,max=255)),cex.lab=1.4,cex.axis=1.3,cex.main=2,space=0.9,font.lab=2)`);
	$R->run(qq`legend("topright",legend=name,pch=15,cex=1.3,bty="n",col=c(rgb(186,85,211,max=255),rgb(65,105,225,max=255),rgb(105,105,105,max=255),rgb(154,205,50,max=255),rgb(0,139,139,max=255),rgb(218,165,32,max=255)))`);
	$R->run(qq`text(bar, number,number,adj=c(100,100),cex=1.3,font=2,pos=3, )`);
	$R->run(qq`dev.off()`);
	$R->stop();
}

if($type =~ /INDEL/i){
	my @len;
	my @indels;
	my @insertions;
	my @deletions;
	my $flag = 1;
	for(my $i = 0; $i<scalar(@indel) and $flag <=20; $i++){
		unless($indel[$i] and $insertion[$i] and $deletion[$i]){
			next;
		}
		$flag++;
		push @len, $i;
		push @indels, $indel[$i]?$indel[$i]:0;
		push @insertions, $insertion[$i]?$insertion[$i]:0;
		push @deletions, $deletion[$i]?$deletion[$i]:0;
	}
	#print "@indels\n";
	#print "@insertions\n";
	#print "@deletions\n";

	#my $dir = dirname $stat;
	my $picture_indel = $dir."/".$sample_name."\.InDel_Distribution.png";
	my $R = Statistics::R->new();
	$R->set('Len', [@len]);
	$R->set('Indel', [@indels]);
	$R->set('Insertion', [@insertions]);
	$R->set('Deletion', [@deletions]);
	$R->run(qq`options(bitmapType='cairo')`);
	$R->run(qq`png("$picture_indel",width=800,heigh=500)`);
	$R->run(qq`par(mar=c(5,5,5,5),mgp=c(3.5,1,0))`);
	$R->run(qq`barplot(rbind(Insertion,Deletion,Indel),beside=T,names.arg=Len,col=c(rgb(255,102,204,max=255),rgb(51,204,204,max=255),rgb(0,102,204,max=255)),ylim=c(0,max(Indel)*1.2),cex.lab=1.5,font.lab=2,cex.axis=1.2,cex.main=2,las=1,xlab="InDel length(bp)",ylab="Number",main="InDel length distribution (All)")`);
	$R->run(qq`legend("right",c("Insertion","Deletion","Indels"),lwd=3,lty=1,bty="n",cex=1.5,col=c(rgb(255,102,204,max=255),rgb(51,204,204,max=255),rgb(0,102,204,max=255)))`);
	$R->run(qq`dev.off()`);
}
