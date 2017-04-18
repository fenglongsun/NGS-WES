#!/bin/bash
#PBS -N lvxy_big_sample.pbs
#PBS -l nodes=1:ppn=4
#PBS -o lvxy_big_sample.log
#PBS -e lvxy_big_sample.err
#PBS -l walltime=1000:00:00

GATK=/opt/NfsDir/BioDir/GenomeAnalysisTK.jar
picard=/opt/NfsDir/BioDir/picard-tools-1.119/picard.jar
samtools=/opt/NfsDir/BioDir/samtools/samtools
igvtools=/opt/NfsDir/BioDir/igvtools_2.3.68/IGVTools/igvtools
script=/opt/NfsDir/BioDir/script/WES_script

ref=/opt/NfsDir/PublicDir/reference/ucsc.hg19.fasta
bait=/opt/NfsDir/PublicDir/AgilentTarget/Bait.interval_list
target=/opt/NfsDir/PublicDir/AgilentTarget/Target.interval_list
EAS_1000g=/opt/NfsDir/PublicDir/1000g/EAS_1000g_Exon.vcf
known1=/opt/NfsDir/PublicDir/1000g/1000G_phase1.indels.hg19.vcf
known2=/opt/NfsDir/PublicDir/1000g/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf
dbsnp=/opt/NfsDir/PublicDir/dbsnp/dbsnp_138.hg19.vcf
dbsnp1=/opt/NfsDir/PublicDir/dbsnp/dbsnp_138.hg19.excluding_sites_after_129.vcf
known3=/opt/NfsDir/NfsDir/PublicDir/1000g/1000G_phase1.snps.high_confidence.hg19.vcf
FASTX=/opt/NfsDir/PublicDir/FASTX
annovar=/opt/NfsDir/PublicDir/annovar

###define a dictory of all sample name
declare -A dic
dic=([1839_1]="1839-1_GATCAG_L008" [1869_2]="1869-2_TAGCTT_L008")
DataPath=/opt/NfsDir/UserDir/lvxy/test/big_sample
for key in $(echo ${!dic[*]})
do

  name=$key
  mkdir  -p $DataPath/$name/output
  OutPath=$DataPath/$name/output
 
  fq1=$DataPath/${dic[$key]}_R1_001.fastq.gz
  fq2=$DataPath/${dic[$key]}_R2_001.fastq.gz


  echo -e Time for quality control of fastq file was started at `date` "\n" >> $OutPath/$name.log

  gzip -dc $fq1  |  $FASTX/bin/fastx_quality_stats -i -  -Q 33 -o $OutPath/$name.R1.stats
  gzip -dc $fq2  |  $FASTX/bin/fastx_quality_stats -i -  -Q 33 -o $OutPath/$name.R2.stats

  awk -F "\t" 'FNR>1{OFS="\t";print $1,$13/($13+$14+$15+$16+$17),$14/($13+$14+$15+$16+$17),$15/($13+$14+$15+$16+$17),$16/($13+$14+$15+$16+$17),$17/($13+$14+$15+$16+$17)}' $OutPath/$name.R1.stats > $OutPath/$name.R1
  Rscript $script/lineplot.R  $OutPath/$name.R1  $OutPath/$name.Nu_distR1.png
  awk -F "\t" 'FNR>1{OFS="\t";print $1,$13/($13+$14+$15+$16+$17),$14/($13+$14+$15+$16+$17),$15/($13+$14+$15+$16+$17),$16/($13+$14+$15+$16+$17),$17/($13+$14+$15+$16+$17)}' $OutPath/$name.R2.stats > $OutPath/$name.R2
  Rscript $script/lineplot.R  $OutPath/$name.R2  $OutPath/$name.Nu_distR2.png

  $FASTX/bin/fastq_quality_boxplot_graph.sh -i $OutPath/$name.R1.stats  -o $OutPath/$name.QC_disR1.png -t $name.R1
  $FASTX/bin/fastq_quality_boxplot_graph.sh -i $OutPath/$name.R2.stats  -o $OutPath/$name.QC_disR2.png -t $name.R2

  echo -e quality control of fastq file was done at `date` "\n" >> $OutPath/$name.log

  echo -e StartTime is `date` "\n" >> $OutPath/$name.log
  bwa mem -t 8  -M -R "@RG\tID:QSY$name\tLB:QSY$name\tSM:$name\tPL:ILLUMINA" $ref $fq1 $fq2 |gzip -3 > $OutPath/$name.align.sam
  echo -e bwa mem done is `date` "\n" >> $OutPath/$name.log

  $samtools view  -@ 8 -bS $OutPath/$name.align.sam -o $OutPath/$name.align.bam
  echo -e samtools view done is `date` "\n" >> $OutPath/$name.log 

  java -jar $picard ReorderSam I= $OutPath/$name.align.bam O= $OutPath/$name.align.reorder.bam R=$ref VALIDATION_STRINGENCY=LENIENT
  echo -e ReorderSamjar done is `date` "\n" >> $OutPath/$name.log

  $samtools sort  $OutPath/$name.align.reorder.bam -o $OutPath/$name.align.reorder.sorted.bam
  echo -e samtools sort done is `date` "\n" >> $OutPath/$name.log

  $samtools index $OutPath/$name.align.reorder.sorted.bam
  echo -e samtools index done is `date` "\n" >> $OutPath/$name.log

  java -jar $picard MarkDuplicates I= $OutPath/$name.align.reorder.sorted.bam O= $OutPath/$name.align.reorder.sorted.markdup.bam ASSUME_SORTED=true METRICS_FILE= $OutPath/$name.align.metrics VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true
  echo -e MarkDuplicatesjar done is `date` "\n" >> $OutPath/$name.log

  java -jar $GATK -T RealignerTargetCreator -R $ref -I $OutPath/$name.align.reorder.sorted.markdup.bam -known $known1 -known $known2 -o $OutPath/$name.realigner.intervals
  echo -e RealignerTargetCreator done is `date` "\n" >> $OutPath/$name.log

  java -jar $GATK -T IndelRealigner -R $ref -I $OutPath/$name.align.reorder.sorted.markdup.bam  -known $known1 -known $known2 -targetIntervals $OutPath/$name.realigner.intervals -o $OutPath/$name.realigner.bam
  echo -e IndelRealigner done is `date` "\n" >> $OutPath/$name.log

  java -jar $GATK -T BaseRecalibrator -R $ref -I $OutPath/$name.realigner.bam  -knownSites $dbsnp -knownSites $known1 -knownSites $known2 -o $OutPath/$name.recal.table 
  echo -e BaseRecalibrator done is `date` "\n" >> $OutPath/$name.log

  java -jar $GATK -T PrintReads -R $ref -I $OutPath/$name.realigner.bam  -BQSR $OutPath/$name.recal.table -o $OutPath/$name.recal.bam
  echo -e PrintReads done is `date` "\n" >> $OutPath/$name.log

  java -jar $GATK -T BaseRecalibrator -R $ref -I $OutPath/$name.realigner.bam  -knownSites $dbsnp -knownSites $known1 -knownSites $known2 -BQSR $OutPath/$name.recal.table -o $OutPath/$name.after_recal.table 
  echo -e BaseRecalibrator done is `date` "\n" >> $OutPath/$name.log

  java -jar $GATK -T AnalyzeCovariates -R $ref -before $OutPath/$name.recal.table -after $OutPath/$name.after_recal.table -plots $OutPath/$name.recal_plots.pdf
  echo -e AnalyzeCovariates done is `date` "\n" >> $OutPath/$name.log

#####verifyBamID
  verifyBamID --vcf $EAS_1000g --bam $OutPath/$name.recal.bam  --out $OutPath/$name.qc1  --verbose -ignoreRG

  awk '{if(NR>1)print "SEQ_ID\t"$1"\nREADS\t"$5"\nAVG_DP\t"$6"\nFREEMIX\t"$7}' $OutPath/$name.qc1.selfSM >$OutPath/$name.mix
#####Picard(OxoG)
#java -jar /opt/NfsDir/Bio_1/lvxy/picard.jar  CollectOxoGMetrics I=$OutPath/$name.recal.bam  O=$OutPath/$name.qc2 R=$ref
#####Picard(SequencingArtifact)
#java -jar /opt/NfsDir/Bio_1/lvxy/picard.jar  CollectSequencingArtifactMetrics I=$OutPath/$name.recal.bam  O=$OutPath/$name.qc3  R=$ref

#####Picard(Hs)
  java -jar $picard  CollectHsMetrics I=$OutPath/$name.recal.bam  O=$OutPath/$name.qc4  R=$ref BAIT_INTERVALS=$bait TARGET_INTERVALS=$target
  perl $script/qc_hs_alig.pl $OutPath/$name.qc4 hs  $name >$OutPath/$name.hs
#####AlignmentSummary
  java -jar $picard  CollectAlignmentSummaryMetrics R=$ref I=$OutPath/$name.recal.bam  O=$OutPath/$name.qc5
  perl $script/qc_hs_alig.pl $OutPath/$name.qc5 alig  $name >$OutPath/$name.alig

  java -jar $picard  CollectQualityYieldMetrics R=$ref I=$OutPath/$name.recal.bam  O=$OutPath/$name.qualigyYield 
  perl $script/qc_hs_alig.pl  $OutPath/$name.qualigyYield qc $name >$OutPath/$name.fq_qc
#####Insertsize
  java -jar $picard  CollectInsertSizeMetrics R=$ref I=$OutPath/$name.recal.bam  O=$OutPath/$name.qc6.txt  H=$OutPath/$name.qc6.pdf M=0.5

  echo -e QC done is `date` "\n" >> $OutPath/$name.log

#java -jar $GATK -T HaplotypeCaller -R $ref -I $OutPath/$name.recal.bam -o $OutPath/$name.g.vcf -ERC GVCF

  java -jar $GATK -T HaplotypeCaller -R $ref -I $OutPath/$name.recal.bam -L $bait -o $OutPath/$name.raw.vcf -stand_call_conf 50 -stand_emit_conf 10.0  -A RMSMappingQuality -A BaseCounts

## harder filering of variant sites
##Extract the SNPs from the call set
  java -jar $GATK  -T SelectVariants -R $ref  -V $OutPath/$name.raw.vcf  -selectType SNP -o $OutPath/$name.raw_snp.vcf

  ## Determine parameters for filtering SNPs 
  ## QD 2.0; FS 60.0; MQ 40.0; MQRankSum -12.5; ReadPosRankSum -8.0; DP 10
  ## Apply the filter to the SNP call set
  java -jar  $GATK -T VariantFiltration -R $ref -V $OutPath/$name.raw_snp.vcf  --filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 ||DP < 10 "  --filterName "my_filter" -o  $OutPath/$name.filtered_snps.vcf

  ## Extract the Indels from the call set
  java -jar $GATK  -T SelectVariants -R $ref  -V $OutPath/$name.raw.vcf  -selectType INDEL -o $OutPath/$name.raw_indel.vcf


  ### Determine parameters for filtering Indels.
  ### QD 2.0; FS 200.0; ReadPosRankSum 20.0; DP 10
  ### Apply the filter to the Indel call set
  java -jar $GATK -T VariantFiltration -R $ref -V  $OutPath/$name.raw_indel.vcf  --filterExpression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || DP < 10"  --filterName "my_filter"  -o $OutPath/$name.filtered_indels.vcf
  grep -v '^#' $OutPath/$name.filtered_indels.vcf >> $OutPath/$name.filtered_snps.vcf
  grep -v  'my_filter' $OutPath/$name.filtered_snps.vcf >$OutPath/$name.final.vcf

  #perl  $script/restrict2bed.pl $bait $OutPath/$name.filter1.vcf > $OutPath/$name.final.vcf
  $igvtools index  $OutPath/$name.final.vcf

  #java -jar $picard  CollectVariantCallingMetrics I=$OutPath/$name.final.vcf O=$OutPath/$name.vcf.dbsnp  DBSNP= $dbsnp TI= $bait 
  #java -jar $picard  CollectVariantCallingMetrics I=$OutPath/$name.final.vcf O=$OutPath/$name.vcf.1000g  DBSNP= $known3  TI= $bait 
  #java -jar $GATK -T VariantEval -R $ref --eval $OutPath/$name.final.vcf   -D $dbsnp1 -noST -ST FunctionalClass -ST Sample  -noEV -EV CompOverlap -EV IndelSummary -EV TiTvVariantEvaluator -EV CountVariants -EV MultiallelicSummary -o $OutPath/$name.eval.report
  echo -e HaplotypeCaller done is `date` "\n" >> $OutPath/$name.log

  #GenotypingGVCFs
  #java -jar $GATK -T GenotypeGVCFs --disable_auto_index_creation_and_locking_when_reading_rods -R $ref -o $OutPath/$name.vcf -D $dbsnp  -L $bait -V $OutPath/$name.g.vcf

  perl $annovar/convert2annovar.pl --includeinfo --format vcf4 --allsample --outfile $OutPath/annovar $OutPath/$name.final.vcf

  perl $annovar/table_annovar.pl --buildver hg19 --thread 3 --remove  --otherinfo --protocol  refGene,cytoBand,avsnp147,clinvar_20160302,cosmic70,1000g2015aug_all,1000g2015aug_afr,1000g2015aug_amr,1000g2015aug_eas,1000g2015aug_eur,1000g2015aug_sas,esp6500siv2_all,esp6500siv2_ea,esp6500siv2_aa,exac03,popfreq_max_20150413,cadd13gt20,gerp++gt2,dbnsfp30a,wgRna,targetScanS,genomicSuperDups,dgvMerged,gwasCatalog,tfbsConsSites,phastConsElements46way -operation g,r,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,r,r,r,r,r,r,r -nastring . $OutPath/annovar.$name.avinput $annovar/humandb --outfile $OutPath/$name.annovar > $OutPath/$name.annovar.log

  head -n 1 $OutPath/$name.annovar.hg19_multianno.txt |awk -F "\t" '{OFS="\t";for (i=1;i<NF;i++) printf("%s\t",$i); print "ref_DP","alt_DP","DP"}'  > $OutPath/$name.annovar.txt
  awk -F "\t" 'FNR>1{OFS="\t";split($90,a,/[:]/);split(a[2],b,/[,]/);for (i=1;i<=80;i++) {printf("%s\t",$i);} {print b[1],b[2],a[3]}}'  $OutPath/$name.annovar.hg19_multianno.txt >> $OutPath/$name.annovar.txt

  perl $script/annovar_statistic.Tp.pl  -a $OutPath/$name.annovar.hg19_multianno.txt -s $name -t SNP
  perl $script/annovar_statistic.Tp.pl  -a $OutPath/$name.annovar.hg19_multianno.txt -s $name -t INDEL
  #paste -d '\t'   $OutPath/$name.SNP_statistic.xls  $OutPath/$name.INDEL_statistic.xls |awk -F "\t" '{print $1"\t"$2"\t"$4}' >$OutPath/$name.vcf_statistic.xls
done

OutPath=$DataPath

echo -e "Sample\t" > head.txt
m=1
for file in `ls $OutPath/*/output/*.SNP_statistic.xls`
do
  mm=${file##*/}
  nn=${mm%%.SNP_statistic.xls}
  sed -i 's/$/&'$nn'\t/' head.txt

  if [ $m == 1 ]
  then
     cp $file  tmp1
     echo $file
     m+=1
  else
     awk -F "\t" '{print $2}' $file > tmp2
     paste -d '\t' tmp1 tmp2|awk -F "\t" '{OFS="\t";print $0}' >> tmp3
     mv tmp3 tmp1
     echo $file
   fi
done

sed -n '1,28p' tmp1 >> head.txt
mv head.txt $OutPath/all_sample.SNP_statistic.xls
rm tmp*

echo -e "Sample\t" > head.txt
m=1
for file in `ls $OutPath/*/output/*.INDEL_statistic.xls`
do
  mm=${file##*/}
  nn=${mm%%.INDEL_statistic.xls}
  sed -i 's/$/&'$nn'\t/' head.txt
  if [ $m == 1 ]
  then
     cp $file  tmp1
     echo $file
     m+=1
  else
     awk -F "\t" '{print $2}' $file > tmp2
     paste -d '\t' tmp1 tmp2|awk -F "\t" '{OFS="\t";print $0}' >> tmp3
     mv tmp3 tmp1
     echo $file
   fi
done

sed -n '1,28p' tmp1 >> head.txt
mv head.txt $OutPath/all_sample.INDEL_statistic.xls
rm tmp*


echo -e "sample\tdownstream\texonic\tintergenic\tintronic\tncRNA_exonic\tncRNA_intronic\tsplicing\tupstream\tUTR3\tUTR5" > $OutPath/summary_figure2.txt;
echo -e "sample\tframeshiftdeletion\tframeshiftinsertion\tnonframeshiftdeletion\tnonframeshiftinsertion\tnonsynonymousSNV\tstopgain\tstoploss\tsynonymous\tunknown" > $OutPath/summary_figure3.txt;

for file in `ls $OutPath/*/output/*annovar.hg19_multianno.txt` 
do
  mm=${file##*/}
  nn=${mm%%.annovar.hg19_multianno.txt}
  echo $nn
  awk -v name=$nn 'BEGIN{sum_downstream=0;sum_exonic=0;sum_intergenic=0;sum_intronic=0;sum_ncRNA_exonic=0;sum_ncRNA_intronic=0;sum_splicing=0;sum_upstream=0;sum_UTR3=0;sum_UTR5=0;sum_other=0}NR>1{OFS="\t";if($6~"downstream"){sum_downstream+=1;}else if($6~"exonic"){sum_exonic+=1;} else if($6~"intergenic"){sum_intergenic+=1;} else if($6~"intronic"){sum_intronic+=1;} else if($6~"ncRNA_exonic"){sum_ncRNA_exonic+=1;} else if($6~"ncRNA_intronic"){sum_ncRNA_intronic+=1;} else if($6~"splicing"){sum_splicing+=1;} else if($6~"upstream"){sum_upstream+=1;} else if($6~"UTR3"){sum_UTR3+=1;} else if($6~"UTR5"){sum_UTR5+=1;} else {sum_other+=1;}}END{print name,sum_downstream,sum_exonic,sum_intergenic,sum_intronic,sum_ncRNA_exonic,sum_ncRNA_intronic,sum_splicing,sum_upstream,sum_UTR3,sum_UTR5}' $OutPath/$nn/output/$nn.annovar.hg19_multianno.txt >> $OutPath/summary_figure2.txt
   awk -v name=$nn 'BEGIN{sum_1=0;sum_2=0;sum_3=0;sum_4=0;sum_5=0;sum_6=0;sum_7=0;sum_8=0;sum_9=0;sum_other=0}NR>1{OFS="\t";if(($9=="frameshift")&&($10=="deletion")){sum_1+=1;}else if(($9=="frameshift")&&($10=="insertion")){sum_2+=1;} else if(($9=="nonframeshift")&&($10=="deletion")){sum_3+=1;} else if(($9=="nonframeshift")&&($10=="insertion")){sum_4+=1;} else if(($9=="nonsynonymous")&&($10=="SNV")){sum_5+=1;} else if($9=="stopgain"){sum_6+=1;} else if($9=="stoploss"){sum_7+=1;} else if(($9=="synonymous")&&($10=="SNV")){sum_8+=1;} else if($9=="unknown"){sum_9+=1;}  else {sum_other+=1;}}END{print name,sum_3,sum_2,sum_3,sum_4,sum_5,sum_6,sum_7,sum_8,sum_9}' $OutPath/$nn/output/$nn.annovar.hg19_multianno.txt >>    $OutPath/summary_figure3.txt
done

Rscript $script/barplot.R  $OutPath/summary_figure2.txt $OutPath/summary_figure2.pdf;
Rscript $script/barplot.R $OutPath/summary_figure3.txt     $OutPath/summary_figure3.pdf


