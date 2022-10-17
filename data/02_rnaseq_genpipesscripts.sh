#!/bin/bash
# Exit immediately on error

set -eu -o pipefail

#-------------------------------------------------------------------------------
# RnaSeqDeNovoAssembly SLURM Job Submission Bash script
# Version: 4.3.1
# Created on: 2022-10-12T13:42:20
# Steps:
#   picard_sam_to_fastq: 0 job... skipping
#   trimmomatic: 4 jobs
#   merge_trimmomatic_stats: 1 job
#   insilico_read_normalization_readsets: 4 jobs
#   insilico_read_normalization_all: 2 jobs
#   trinity: 2 jobs
#   exonerate_fastasplit: 1 job
#   blastx_trinity_uniprot: 20 jobs
#   blastx_trinity_uniprot_merge: 2 jobs
#   transdecoder: 1 job
#   hmmer: 1 job
#   rnammer_transcriptome: 1 job
#   blastp_transdecoder_uniprot: 1 job
#   signalp: 1 job
#   tmhmm: 1 job
#   trinotate: 2 jobs
#   align_and_estimate_abundance_prep_reference: 1 job
#   align_and_estimate_abundance: 7 jobs
#   gq_seq_utils_exploratory_analysis_rnaseq_denovo: 2 jobs
#   differential_expression: 3 jobs
#   filter_annotated_components: 2 jobs
#   gq_seq_utils_exploratory_analysis_rnaseq_denovo_filtered: 3 jobs
#   differential_expression_filtered: 4 jobs
#   TOTAL: 66 jobs
#-------------------------------------------------------------------------------

OUTPUT_DIR=/lustre04/scratch/senthil/rnaseq
JOB_OUTPUT_DIR=$OUTPUT_DIR/job_output
TIMESTAMP=`date +%FT%H.%M.%S`
JOB_LIST=$JOB_OUTPUT_DIR/RnaSeqDeNovoAssembly_job_list_$TIMESTAMP
export CONFIG_FILES="/cvmfs/soft.mugqic/root/software/genpipes/genpipes-4.3.1/pipelines/rnaseq_denovo_assembly/rnaseq_denovo_assembly.base.ini,/cvmfs/soft.mugqic/root/software/genpipes/genpipes-4.3.1/pipelines/common_ini/beluga.ini"
mkdir -p $OUTPUT_DIR
cd $OUTPUT_DIR

#-------------------------------------------------------------------------------
# STEP: trimmomatic
#-------------------------------------------------------------------------------
STEP=trimmomatic
mkdir -p $JOB_OUTPUT_DIR/$STEP

#-------------------------------------------------------------------------------
# JOB: trimmomatic_1_JOB_ID: trimmomatic.GM12878_Rep1
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.GM12878_Rep1
JOB_DEPENDENCIES=
JOB_DONE=job_output/trimmomatic/trimmomatic.GM12878_Rep1.6e020bb72580524636df16d3d09ccd96.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'trimmomatic.GM12878_Rep1.6e020bb72580524636df16d3d09ccd96.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.7.0_60 mugqic/trimmomatic/0.36 && \
mkdir -p trim/GM12878_Rep1 && \
`cat > trim/GM12878_Rep1/GM12878_Rep1.trim.adapters.fa << END
>Prefix/1
ACACTCTTTCCCTACACGACGCTCTTCCGATCT
>Prefix/2
TGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
END
` && \
java -XX:ParallelGCThreads=4 -Xmx10000M -jar $TRIMMOMATIC_JAR PE \
  -threads 10 \
  -phred33 \
  /cvmfs/soft.mugqic/CentOS6/testdata/rnaseq/raw_data/rnaseq_GM12878_chr19_Rep1_R1.fastq.gz \
  /cvmfs/soft.mugqic/CentOS6/testdata/rnaseq/raw_data/rnaseq_GM12878_chr19_Rep1_R2.fastq.gz \
  trim/GM12878_Rep1/GM12878_Rep1.trim.pair1.fastq.gz \
  trim/GM12878_Rep1/GM12878_Rep1.trim.single1.fastq.gz \
  trim/GM12878_Rep1/GM12878_Rep1.trim.pair2.fastq.gz \
  trim/GM12878_Rep1/GM12878_Rep1.trim.single2.fastq.gz \
  ILLUMINACLIP:trim/GM12878_Rep1/GM12878_Rep1.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/GM12878_Rep1/GM12878_Rep1.trim.log
trimmomatic.GM12878_Rep1.6e020bb72580524636df16d3d09ccd96.mugqic.done
chmod 755 $COMMAND
trimmomatic_1_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&     $COMMAND 
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE


if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem 10G -c 10 -N 1    | grep "[0-9]" | cut -d\  -f4)
echo "$trimmomatic_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$trimmomatic_1_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: trimmomatic_2_JOB_ID: trimmomatic.GM12878_Rep2
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.GM12878_Rep2
JOB_DEPENDENCIES=
JOB_DONE=job_output/trimmomatic/trimmomatic.GM12878_Rep2.3165861387368aabc58ec2c887505d73.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'trimmomatic.GM12878_Rep2.3165861387368aabc58ec2c887505d73.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.7.0_60 mugqic/trimmomatic/0.36 && \
mkdir -p trim/GM12878_Rep2 && \
`cat > trim/GM12878_Rep2/GM12878_Rep2.trim.adapters.fa << END
>Prefix/1
ACACTCTTTCCCTACACGACGCTCTTCCGATCT
>Prefix/2
TGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
END
` && \
java -XX:ParallelGCThreads=4 -Xmx10000M -jar $TRIMMOMATIC_JAR PE \
  -threads 10 \
  -phred33 \
  /cvmfs/soft.mugqic/CentOS6/testdata/rnaseq/raw_data/rnaseq_GM12878_chr19_Rep2_R1.fastq.gz \
  /cvmfs/soft.mugqic/CentOS6/testdata/rnaseq/raw_data/rnaseq_GM12878_chr19_Rep2_R2.fastq.gz \
  trim/GM12878_Rep2/GM12878_Rep2.trim.pair1.fastq.gz \
  trim/GM12878_Rep2/GM12878_Rep2.trim.single1.fastq.gz \
  trim/GM12878_Rep2/GM12878_Rep2.trim.pair2.fastq.gz \
  trim/GM12878_Rep2/GM12878_Rep2.trim.single2.fastq.gz \
  ILLUMINACLIP:trim/GM12878_Rep2/GM12878_Rep2.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/GM12878_Rep2/GM12878_Rep2.trim.log
trimmomatic.GM12878_Rep2.3165861387368aabc58ec2c887505d73.mugqic.done
chmod 755 $COMMAND
trimmomatic_2_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&     $COMMAND 
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE


if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem 10G -c 10 -N 1    | grep "[0-9]" | cut -d\  -f4)
echo "$trimmomatic_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$trimmomatic_2_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: trimmomatic_3_JOB_ID: trimmomatic.H1ESC_Rep1
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.H1ESC_Rep1
JOB_DEPENDENCIES=
JOB_DONE=job_output/trimmomatic/trimmomatic.H1ESC_Rep1.375487d7453fe48924419ead231cb595.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'trimmomatic.H1ESC_Rep1.375487d7453fe48924419ead231cb595.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.7.0_60 mugqic/trimmomatic/0.36 && \
mkdir -p trim/H1ESC_Rep1 && \
`cat > trim/H1ESC_Rep1/H1ESC_Rep1.trim.adapters.fa << END
>Prefix/1
ACACTCTTTCCCTACACGACGCTCTTCCGATCT
>Prefix/2
TGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
END
` && \
java -XX:ParallelGCThreads=4 -Xmx10000M -jar $TRIMMOMATIC_JAR PE \
  -threads 10 \
  -phred33 \
  /cvmfs/soft.mugqic/CentOS6/testdata/rnaseq/raw_data/rnaseq_H1ESC_chr19_Rep1_R1.fastq.gz \
  /cvmfs/soft.mugqic/CentOS6/testdata/rnaseq/raw_data/rnaseq_H1ESC_chr19_Rep1_R2.fastq.gz \
  trim/H1ESC_Rep1/H1ESC_Rep1.trim.pair1.fastq.gz \
  trim/H1ESC_Rep1/H1ESC_Rep1.trim.single1.fastq.gz \
  trim/H1ESC_Rep1/H1ESC_Rep1.trim.pair2.fastq.gz \
  trim/H1ESC_Rep1/H1ESC_Rep1.trim.single2.fastq.gz \
  ILLUMINACLIP:trim/H1ESC_Rep1/H1ESC_Rep1.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/H1ESC_Rep1/H1ESC_Rep1.trim.log
trimmomatic.H1ESC_Rep1.375487d7453fe48924419ead231cb595.mugqic.done
chmod 755 $COMMAND
trimmomatic_3_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&     $COMMAND 
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE


if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem 10G -c 10 -N 1    | grep "[0-9]" | cut -d\  -f4)
echo "$trimmomatic_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$trimmomatic_3_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: trimmomatic_4_JOB_ID: trimmomatic.H1ESC_Rep2
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.H1ESC_Rep2
JOB_DEPENDENCIES=
JOB_DONE=job_output/trimmomatic/trimmomatic.H1ESC_Rep2.2769354918ff10ba289cbb839e944c95.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'trimmomatic.H1ESC_Rep2.2769354918ff10ba289cbb839e944c95.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.7.0_60 mugqic/trimmomatic/0.36 && \
mkdir -p trim/H1ESC_Rep2 && \
`cat > trim/H1ESC_Rep2/H1ESC_Rep2.trim.adapters.fa << END
>Prefix/1
ACACTCTTTCCCTACACGACGCTCTTCCGATCT
>Prefix/2
TGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
END
` && \
java -XX:ParallelGCThreads=4 -Xmx10000M -jar $TRIMMOMATIC_JAR PE \
  -threads 10 \
  -phred33 \
  /cvmfs/soft.mugqic/CentOS6/testdata/rnaseq/raw_data/rnaseq_H1ESC_chr19_Rep2_R1.fastq.gz \
  /cvmfs/soft.mugqic/CentOS6/testdata/rnaseq/raw_data/rnaseq_H1ESC_chr19_Rep2_R2.fastq.gz \
  trim/H1ESC_Rep2/H1ESC_Rep2.trim.pair1.fastq.gz \
  trim/H1ESC_Rep2/H1ESC_Rep2.trim.single1.fastq.gz \
  trim/H1ESC_Rep2/H1ESC_Rep2.trim.pair2.fastq.gz \
  trim/H1ESC_Rep2/H1ESC_Rep2.trim.single2.fastq.gz \
  ILLUMINACLIP:trim/H1ESC_Rep2/H1ESC_Rep2.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/H1ESC_Rep2/H1ESC_Rep2.trim.log
trimmomatic.H1ESC_Rep2.2769354918ff10ba289cbb839e944c95.mugqic.done
chmod 755 $COMMAND
trimmomatic_4_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&     $COMMAND 
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE


if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem 10G -c 10 -N 1    | grep "[0-9]" | cut -d\  -f4)
echo "$trimmomatic_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$trimmomatic_4_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# STEP: merge_trimmomatic_stats
#-------------------------------------------------------------------------------
STEP=merge_trimmomatic_stats
mkdir -p $JOB_OUTPUT_DIR/$STEP

#-------------------------------------------------------------------------------
# JOB: merge_trimmomatic_stats_1_JOB_ID: merge_trimmomatic_stats
#-------------------------------------------------------------------------------
JOB_NAME=merge_trimmomatic_stats
JOB_DEPENDENCIES=$trimmomatic_1_JOB_ID:$trimmomatic_2_JOB_ID:$trimmomatic_3_JOB_ID:$trimmomatic_4_JOB_ID
JOB_DONE=job_output/merge_trimmomatic_stats/merge_trimmomatic_stats.2c3d567ecea3ffdffcc18075f5b6fc62.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'merge_trimmomatic_stats.2c3d567ecea3ffdffcc18075f5b6fc62.mugqic.done' > $COMMAND
module purge && \
module load mugqic/perl/5.22.1 mugqic/pandoc/2.16.1 && \
mkdir -p metrics && \
echo 'Sample	Readset	Raw Paired Reads #	Surviving Paired Reads #	Surviving Paired Reads %' > metrics/trimReadsetTable.tsv && \
grep ^Input trim/GM12878_Rep1/GM12878_Rep1.trim.log | \
perl -pe 's/^Input Read Pairs: (\d+).*Both Surviving: (\d+).*Forward Only Surviving: (\d+).*$/GM12878_Rep1	GM12878_Rep1	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/GM12878_Rep2/GM12878_Rep2.trim.log | \
perl -pe 's/^Input Read Pairs: (\d+).*Both Surviving: (\d+).*Forward Only Surviving: (\d+).*$/GM12878_Rep2	GM12878_Rep2	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/H1ESC_Rep1/H1ESC_Rep1.trim.log | \
perl -pe 's/^Input Read Pairs: (\d+).*Both Surviving: (\d+).*Forward Only Surviving: (\d+).*$/H1ESC_Rep1	H1ESC_Rep1	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/H1ESC_Rep2/H1ESC_Rep2.trim.log | \
perl -pe 's/^Input Read Pairs: (\d+).*Both Surviving: (\d+).*Forward Only Surviving: (\d+).*$/H1ESC_Rep2	H1ESC_Rep2	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
cut -f1,3- metrics/trimReadsetTable.tsv | awk -F"	" '{OFS="	"; if (NR==1) {if ($2=="Raw Paired Reads #") {paired=1};print "Sample", "Raw Reads #", "Surviving Reads #", "Surviving %"} else {if (paired) {$2=$2*2; $3=$3*2}; raw[$1]+=$2; surviving[$1]+=$3}}END{for (sample in raw){print sample, raw[sample], surviving[sample], surviving[sample] / raw[sample] * 100}}' \
  > metrics/trimSampleTable.tsv && \
mkdir -p report && \
cp metrics/trimReadsetTable.tsv metrics/trimSampleTable.tsv report/ && \
trim_readset_table_md=`LC_NUMERIC=en_CA awk -F "	" '{OFS="|"; if (NR == 1) {$1 = $1; print $0; print "-----|-----|-----:|-----:|-----:"} else {print $1, $2, sprintf("%\47d", $3), sprintf("%\47d", $4), sprintf("%.1f", $5)}}' metrics/trimReadsetTable.tsv` && \
pandoc \
  /cvmfs/soft.mugqic/root/software/genpipes/genpipes-4.3.1/bfx/report/Illumina.merge_trimmomatic_stats.md \
  --template /cvmfs/soft.mugqic/root/software/genpipes/genpipes-4.3.1/bfx/report/Illumina.merge_trimmomatic_stats.md \
  --variable trailing_min_quality=30 \
  --variable min_length=50 \
  --variable read_type=Paired \
  --variable trim_readset_table="$trim_readset_table_md" \
  --to markdown \
  > report/Illumina.merge_trimmomatic_stats.md
merge_trimmomatic_stats.2c3d567ecea3ffdffcc18075f5b6fc62.mugqic.done
chmod 755 $COMMAND
merge_trimmomatic_stats_1_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&     $COMMAND 
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE


if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem-per-cpu 4700M -c 1 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$merge_trimmomatic_stats_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$merge_trimmomatic_stats_1_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# STEP: insilico_read_normalization_readsets
#-------------------------------------------------------------------------------
STEP=insilico_read_normalization_readsets
mkdir -p $JOB_OUTPUT_DIR/$STEP

#-------------------------------------------------------------------------------
# JOB: insilico_read_normalization_readsets_1_JOB_ID: insilico_read_normalization_readsets.GM12878_Rep1
#-------------------------------------------------------------------------------
JOB_NAME=insilico_read_normalization_readsets.GM12878_Rep1
JOB_DEPENDENCIES=$trimmomatic_1_JOB_ID
JOB_DONE=job_output/insilico_read_normalization_readsets/insilico_read_normalization_readsets.GM12878_Rep1.6a975696f476fbffe9bb2da8bc1badc7.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'insilico_read_normalization_readsets.GM12878_Rep1.6a975696f476fbffe9bb2da8bc1badc7.mugqic.done' > $COMMAND
module purge && \
module load mugqic/perl/5.22.1 mugqic/trinity/2.0.4_patch && \
mkdir -p insilico_read_normalization/GM12878_Rep1 && \
insilico_read_normalization.pl --pairs_together --SS_lib_type RF --PARALLEL_STATS --KMER_SIZE 25 --max_pct_stdev 100 \
  --seqType fq \
  --JM 180G \
  --max_cov 30 \
  --left trim/GM12878_Rep1/GM12878_Rep1.trim.pair1.fastq.gz \
  --right trim/GM12878_Rep1/GM12878_Rep1.trim.pair2.fastq.gz \
  --output insilico_read_normalization/GM12878_Rep1 \
  --CPU 40 && \
wc -l insilico_read_normalization/GM12878_Rep1/left.norm.fq | awk '{print "# normalized paired reads	"$1 / 4}' > insilico_read_normalization/GM12878_Rep1/normalization.stats.tsv
insilico_read_normalization_readsets.GM12878_Rep1.6a975696f476fbffe9bb2da8bc1badc7.mugqic.done
chmod 755 $COMMAND
insilico_read_normalization_readsets_1_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&     $COMMAND 
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE


if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem 187G -c 40 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$insilico_read_normalization_readsets_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$insilico_read_normalization_readsets_1_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: insilico_read_normalization_readsets_2_JOB_ID: insilico_read_normalization_readsets.GM12878_Rep2
#-------------------------------------------------------------------------------
JOB_NAME=insilico_read_normalization_readsets.GM12878_Rep2
JOB_DEPENDENCIES=$trimmomatic_2_JOB_ID
JOB_DONE=job_output/insilico_read_normalization_readsets/insilico_read_normalization_readsets.GM12878_Rep2.6a9710a6bb1d335cd32892a59fc7ea10.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'insilico_read_normalization_readsets.GM12878_Rep2.6a9710a6bb1d335cd32892a59fc7ea10.mugqic.done' > $COMMAND
module purge && \
module load mugqic/perl/5.22.1 mugqic/trinity/2.0.4_patch && \
mkdir -p insilico_read_normalization/GM12878_Rep2 && \
insilico_read_normalization.pl --pairs_together --SS_lib_type RF --PARALLEL_STATS --KMER_SIZE 25 --max_pct_stdev 100 \
  --seqType fq \
  --JM 180G \
  --max_cov 30 \
  --left trim/GM12878_Rep2/GM12878_Rep2.trim.pair1.fastq.gz \
  --right trim/GM12878_Rep2/GM12878_Rep2.trim.pair2.fastq.gz \
  --output insilico_read_normalization/GM12878_Rep2 \
  --CPU 40 && \
wc -l insilico_read_normalization/GM12878_Rep2/left.norm.fq | awk '{print "# normalized paired reads	"$1 / 4}' > insilico_read_normalization/GM12878_Rep2/normalization.stats.tsv
insilico_read_normalization_readsets.GM12878_Rep2.6a9710a6bb1d335cd32892a59fc7ea10.mugqic.done
chmod 755 $COMMAND
insilico_read_normalization_readsets_2_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&     $COMMAND 
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE


if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem 187G -c 40 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$insilico_read_normalization_readsets_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$insilico_read_normalization_readsets_2_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: insilico_read_normalization_readsets_3_JOB_ID: insilico_read_normalization_readsets.H1ESC_Rep1
#-------------------------------------------------------------------------------
JOB_NAME=insilico_read_normalization_readsets.H1ESC_Rep1
JOB_DEPENDENCIES=$trimmomatic_3_JOB_ID
JOB_DONE=job_output/insilico_read_normalization_readsets/insilico_read_normalization_readsets.H1ESC_Rep1.36e1988b9a16569527b361358b7e1915.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'insilico_read_normalization_readsets.H1ESC_Rep1.36e1988b9a16569527b361358b7e1915.mugqic.done' > $COMMAND
module purge && \
module load mugqic/perl/5.22.1 mugqic/trinity/2.0.4_patch && \
mkdir -p insilico_read_normalization/H1ESC_Rep1 && \
insilico_read_normalization.pl --pairs_together --SS_lib_type RF --PARALLEL_STATS --KMER_SIZE 25 --max_pct_stdev 100 \
  --seqType fq \
  --JM 180G \
  --max_cov 30 \
  --left trim/H1ESC_Rep1/H1ESC_Rep1.trim.pair1.fastq.gz \
  --right trim/H1ESC_Rep1/H1ESC_Rep1.trim.pair2.fastq.gz \
  --output insilico_read_normalization/H1ESC_Rep1 \
  --CPU 40 && \
wc -l insilico_read_normalization/H1ESC_Rep1/left.norm.fq | awk '{print "# normalized paired reads	"$1 / 4}' > insilico_read_normalization/H1ESC_Rep1/normalization.stats.tsv
insilico_read_normalization_readsets.H1ESC_Rep1.36e1988b9a16569527b361358b7e1915.mugqic.done
chmod 755 $COMMAND
insilico_read_normalization_readsets_3_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&     $COMMAND 
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE


if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem 187G -c 40 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$insilico_read_normalization_readsets_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$insilico_read_normalization_readsets_3_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: insilico_read_normalization_readsets_4_JOB_ID: insilico_read_normalization_readsets.H1ESC_Rep2
#-------------------------------------------------------------------------------
JOB_NAME=insilico_read_normalization_readsets.H1ESC_Rep2
JOB_DEPENDENCIES=$trimmomatic_4_JOB_ID
JOB_DONE=job_output/insilico_read_normalization_readsets/insilico_read_normalization_readsets.H1ESC_Rep2.ef555b414747548b847f8d38583680a8.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'insilico_read_normalization_readsets.H1ESC_Rep2.ef555b414747548b847f8d38583680a8.mugqic.done' > $COMMAND
module purge && \
module load mugqic/perl/5.22.1 mugqic/trinity/2.0.4_patch && \
mkdir -p insilico_read_normalization/H1ESC_Rep2 && \
insilico_read_normalization.pl --pairs_together --SS_lib_type RF --PARALLEL_STATS --KMER_SIZE 25 --max_pct_stdev 100 \
  --seqType fq \
  --JM 180G \
  --max_cov 30 \
  --left trim/H1ESC_Rep2/H1ESC_Rep2.trim.pair1.fastq.gz \
  --right trim/H1ESC_Rep2/H1ESC_Rep2.trim.pair2.fastq.gz \
  --output insilico_read_normalization/H1ESC_Rep2 \
  --CPU 40 && \
wc -l insilico_read_normalization/H1ESC_Rep2/left.norm.fq | awk '{print "# normalized paired reads	"$1 / 4}' > insilico_read_normalization/H1ESC_Rep2/normalization.stats.tsv
insilico_read_normalization_readsets.H1ESC_Rep2.ef555b414747548b847f8d38583680a8.mugqic.done
chmod 755 $COMMAND
insilico_read_normalization_readsets_4_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&     $COMMAND 
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE


if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem 187G -c 40 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$insilico_read_normalization_readsets_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$insilico_read_normalization_readsets_4_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# STEP: insilico_read_normalization_all
#-------------------------------------------------------------------------------
STEP=insilico_read_normalization_all
mkdir -p $JOB_OUTPUT_DIR/$STEP

#-------------------------------------------------------------------------------
# JOB: insilico_read_normalization_all_1_JOB_ID: insilico_read_normalization_all
#-------------------------------------------------------------------------------
JOB_NAME=insilico_read_normalization_all
JOB_DEPENDENCIES=$insilico_read_normalization_readsets_1_JOB_ID:$insilico_read_normalization_readsets_2_JOB_ID:$insilico_read_normalization_readsets_3_JOB_ID:$insilico_read_normalization_readsets_4_JOB_ID
JOB_DONE=job_output/insilico_read_normalization_all/insilico_read_normalization_all.fd423ce965c2ebb0c3be6b28de4b4bdf.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'insilico_read_normalization_all.fd423ce965c2ebb0c3be6b28de4b4bdf.mugqic.done' > $COMMAND
module purge && \
module load mugqic/perl/5.22.1 mugqic/trinity/2.0.4_patch && \
mkdir -p insilico_read_normalization/all && \
insilico_read_normalization.pl --pairs_together --SS_lib_type RF --PARALLEL_STATS --KMER_SIZE 25 --max_pct_stdev 100 \
  --seqType fq \
  --JM 180G \
  --max_cov 30 \
  --left insilico_read_normalization/GM12878_Rep1/left.norm.fq \
  --left insilico_read_normalization/GM12878_Rep2/left.norm.fq \
  --left insilico_read_normalization/H1ESC_Rep1/left.norm.fq \
  --left insilico_read_normalization/H1ESC_Rep2/left.norm.fq \
  --right insilico_read_normalization/GM12878_Rep1/right.norm.fq \
  --right insilico_read_normalization/GM12878_Rep2/right.norm.fq \
  --right insilico_read_normalization/H1ESC_Rep1/right.norm.fq \
  --right insilico_read_normalization/H1ESC_Rep2/right.norm.fq \
  --output insilico_read_normalization/all \
  --CPU 40 && \
wc -l insilico_read_normalization/all/left.norm.fq | awk '{print "# normalized paired reads	"$1 / 4}' > insilico_read_normalization/all/normalization.stats.tsv
insilico_read_normalization_all.fd423ce965c2ebb0c3be6b28de4b4bdf.mugqic.done
chmod 755 $COMMAND
insilico_read_normalization_all_1_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&     $COMMAND 
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE


if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem 187G -c 40 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$insilico_read_normalization_all_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$insilico_read_normalization_all_1_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: insilico_read_normalization_all_2_JOB_ID: insilico_read_normalization_all_report
#-------------------------------------------------------------------------------
JOB_NAME=insilico_read_normalization_all_report
JOB_DEPENDENCIES=$merge_trimmomatic_stats_1_JOB_ID:$insilico_read_normalization_all_1_JOB_ID
JOB_DONE=job_output/insilico_read_normalization_all/insilico_read_normalization_all_report.23abe30984267ba58d41b2490202afd1.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'insilico_read_normalization_all_report.23abe30984267ba58d41b2490202afd1.mugqic.done' > $COMMAND
module purge && \
module load mugqic/pandoc/2.16.1 && \
mkdir -p report && \
sum_norm=`cut -f2 insilico_read_normalization/all/normalization.stats.tsv` && \
normalization_table=`sed '1d' report/trimReadsetTable.tsv | LC_NUMERIC=en_CA awk -v sum_norm=$sum_norm '{sum_trim+=$4}END{print sprintf("%\47d", sum_trim)"|"sprintf("%\47d", sum_norm)"|"sprintf("%.2f", sum_norm / sum_trim * 100)}'` && \
pandoc --to=markdown \
  --template /cvmfs/soft.mugqic/root/software/genpipes/genpipes-4.3.1/bfx/report/RnaSeqDeNovoAssembly.insilico_read_normalization_all.md \
  --variable read_type="Paired" \
  --variable normalization_table="$normalization_table" \
  /cvmfs/soft.mugqic/root/software/genpipes/genpipes-4.3.1/bfx/report/RnaSeqDeNovoAssembly.insilico_read_normalization_all.md \
  > report/RnaSeqDeNovoAssembly.insilico_read_normalization_all.md
insilico_read_normalization_all_report.23abe30984267ba58d41b2490202afd1.mugqic.done
chmod 755 $COMMAND
insilico_read_normalization_all_2_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&     $COMMAND 
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE


if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem-per-cpu 4700M -c 1 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$insilico_read_normalization_all_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$insilico_read_normalization_all_2_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# STEP: trinity
#-------------------------------------------------------------------------------
STEP=trinity
mkdir -p $JOB_OUTPUT_DIR/$STEP

#-------------------------------------------------------------------------------
# JOB: trinity_1_JOB_ID: trinity
#-------------------------------------------------------------------------------
JOB_NAME=trinity
JOB_DEPENDENCIES=$insilico_read_normalization_all_1_JOB_ID
JOB_DONE=job_output/trinity/trinity.62550400241a7e75d4b6c1eecb0fa246.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'trinity.62550400241a7e75d4b6c1eecb0fa246.mugqic.done' > $COMMAND
module purge && \
module load mugqic/perl/5.22.1 mugqic/java/openjdk-jdk1.7.0_60 mugqic/trinity/2.0.4_patch mugqic/bowtie/1.0.0 mugqic/samtools/1.12 mugqic/R_Bioconductor/3.6.0_3.9 mugqic/mugqic_R_packages/1.0.6 && \
Trinity --seqType fq --SS_lib_type RF --min_contig_length 200 --min_kmer_cov 2 --bflyGCThreads 1 --bflyHeapSpaceMax 10G --bflyCPU 8 --no_version_check \
  --max_memory 180G \
  --CPU 40 \
  --left insilico_read_normalization/all/left.norm.fq \
  --right insilico_read_normalization/all/right.norm.fq \
  --output trinity_out_dir && \
zip -j trinity_out_dir/Trinity.fasta.zip trinity_out_dir/Trinity.fasta && \
Rscript -e 'library(gqSeqUtils); dnaFastaStats(filename = "trinity_out_dir/Trinity.fasta", type = "trinity", output.prefix = "trinity_out_dir/Trinity.stats")'
trinity.62550400241a7e75d4b6c1eecb0fa246.mugqic.done
chmod 755 $COMMAND
trinity_1_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&     $COMMAND 
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE


if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem 187G -c 40 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$trinity_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$trinity_1_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: trinity_2_JOB_ID: trinity_report
#-------------------------------------------------------------------------------
JOB_NAME=trinity_report
JOB_DEPENDENCIES=$trinity_1_JOB_ID
JOB_DONE=job_output/trinity/trinity_report.39d5e60eedaafc2e5e4b4146916da129.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'trinity_report.39d5e60eedaafc2e5e4b4146916da129.mugqic.done' > $COMMAND
module purge && \
module load mugqic/pandoc/2.16.1 && \
mkdir -p report && \
cp trinity_out_dir/Trinity.fasta.zip trinity_out_dir/Trinity.stats.csv trinity_out_dir/Trinity.stats.jpg trinity_out_dir/Trinity.stats.pdf report/ && \
assembly_table=`sed '1d' trinity_out_dir/Trinity.stats.csv | perl -pe 's/^"([^"]*)",/\1	/g' | grep -P "^(Nb. Transcripts|Nb. Components|Total Transcripts Length|Min. Transcript Length|Median Transcript Length|Mean Transcript Length|Max. Transcript Length|N50)" | LC_NUMERIC=en_CA awk -F"	" '{print $1"|"sprintf("%\47d", $2)}'` && \
pandoc --to=markdown \
  --template /cvmfs/soft.mugqic/root/software/genpipes/genpipes-4.3.1/bfx/report/RnaSeqDeNovoAssembly.trinity.md \
  --variable assembly_table="$assembly_table" \
  /cvmfs/soft.mugqic/root/software/genpipes/genpipes-4.3.1/bfx/report/RnaSeqDeNovoAssembly.trinity.md \
  > report/RnaSeqDeNovoAssembly.trinity.md
trinity_report.39d5e60eedaafc2e5e4b4146916da129.mugqic.done
chmod 755 $COMMAND
trinity_2_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&     $COMMAND 
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE


if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem-per-cpu 4700M -c 1 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$trinity_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$trinity_2_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# STEP: exonerate_fastasplit
#-------------------------------------------------------------------------------
STEP=exonerate_fastasplit
mkdir -p $JOB_OUTPUT_DIR/$STEP

#-------------------------------------------------------------------------------
# JOB: exonerate_fastasplit_1_JOB_ID: exonerate_fastasplit.Trinity.fasta
#-------------------------------------------------------------------------------
JOB_NAME=exonerate_fastasplit.Trinity.fasta
JOB_DEPENDENCIES=$trinity_1_JOB_ID
JOB_DONE=job_output/exonerate_fastasplit/exonerate_fastasplit.Trinity.fasta.c92b42e79873a8371f77a03cdec520e2.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'exonerate_fastasplit.Trinity.fasta.c92b42e79873a8371f77a03cdec520e2.mugqic.done' > $COMMAND
module purge && \
module load mugqic/exonerate/2.2.0 && \
rm -rf trinity_out_dir/Trinity.fasta_chunks && \
mkdir -p trinity_out_dir/Trinity.fasta_chunks && \
awk '{ print $1 }' trinity_out_dir/Trinity.fasta  > trinity_out_dir/Trinity.fa && \
fastasplit -f trinity_out_dir/Trinity.fa -o trinity_out_dir/Trinity.fasta_chunks -c 20
exonerate_fastasplit.Trinity.fasta.c92b42e79873a8371f77a03cdec520e2.mugqic.done
chmod 755 $COMMAND
exonerate_fastasplit_1_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&     $COMMAND 
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE


if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem-per-cpu 4700M -c 1 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$exonerate_fastasplit_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$exonerate_fastasplit_1_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# STEP: blastx_trinity_uniprot
#-------------------------------------------------------------------------------
STEP=blastx_trinity_uniprot
mkdir -p $JOB_OUTPUT_DIR/$STEP

#-------------------------------------------------------------------------------
# JOB: blastx_trinity_uniprot_1_JOB_ID: blastx_trinity_uniprot.uniprot_sprot.trinotate_v2.0.pep.chunk_0000000
#-------------------------------------------------------------------------------
JOB_NAME=blastx_trinity_uniprot.uniprot_sprot.trinotate_v2.0.pep.chunk_0000000
JOB_DEPENDENCIES=$exonerate_fastasplit_1_JOB_ID
JOB_DONE=job_output/blastx_trinity_uniprot/blastx_trinity_uniprot.uniprot_sprot.trinotate_v2.0.pep.chunk_0000000.7db0b0a8fe98ef8089dc88a1614fefa6.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'blastx_trinity_uniprot.uniprot_sprot.trinotate_v2.0.pep.chunk_0000000.7db0b0a8fe98ef8089dc88a1614fefa6.mugqic.done' > $COMMAND
module purge && \
module load mugqic/perl/5.22.1 mugqic/mugqic_tools/2.8.2 mugqic/blast/2.3.0+ && \
mkdir -p blast && \
ln -s -f ../trinity_out_dir/Trinity.fasta_chunks/Trinity.fa_chunk_0000000 blast/query_Trinity_uniprot_sprot.trinotate_v2.0.pep_chunk_0000000.tsv && \
parallelBlast.pl \
-file blast/query_Trinity_uniprot_sprot.trinotate_v2.0.pep_chunk_0000000.tsv \
--OUT blast/blastx_Trinity_uniprot_sprot.trinotate_v2.0.pep_chunk_0000000.tsv \
-n 40 \
--BLAST "blastx -db /cvmfs/soft.mugqic/CentOS6/genomes/blast_db/uniprot_sprot.trinotate_v2.0.pep -max_target_seqs 1 -outfmt '6 std stitle'" 
blastx_trinity_uniprot.uniprot_sprot.trinotate_v2.0.pep.chunk_0000000.7db0b0a8fe98ef8089dc88a1614fefa6.mugqic.done
chmod 755 $COMMAND
blastx_trinity_uniprot_1_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&     $COMMAND 
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE


if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem-per-cpu 4700M -c 40 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$blastx_trinity_uniprot_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$blastx_trinity_uniprot_1_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: blastx_trinity_uniprot_2_JOB_ID: blastx_trinity_uniprot.uniprot_sprot.trinotate_v2.0.pep.chunk_0000001
#-------------------------------------------------------------------------------
JOB_NAME=blastx_trinity_uniprot.uniprot_sprot.trinotate_v2.0.pep.chunk_0000001
JOB_DEPENDENCIES=$exonerate_fastasplit_1_JOB_ID
JOB_DONE=job_output/blastx_trinity_uniprot/blastx_trinity_uniprot.uniprot_sprot.trinotate_v2.0.pep.chunk_0000001.aa0c502502949b69d484967362fdb1a8.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'blastx_trinity_uniprot.uniprot_sprot.trinotate_v2.0.pep.chunk_0000001.aa0c502502949b69d484967362fdb1a8.mugqic.done' > $COMMAND
module purge && \
module load mugqic/perl/5.22.1 mugqic/mugqic_tools/2.8.2 mugqic/blast/2.3.0+ && \
mkdir -p blast && \
ln -s -f ../trinity_out_dir/Trinity.fasta_chunks/Trinity.fa_chunk_0000001 blast/query_Trinity_uniprot_sprot.trinotate_v2.0.pep_chunk_0000001.tsv && \
parallelBlast.pl \
-file blast/query_Trinity_uniprot_sprot.trinotate_v2.0.pep_chunk_0000001.tsv \
--OUT blast/blastx_Trinity_uniprot_sprot.trinotate_v2.0.pep_chunk_0000001.tsv \
-n 40 \
--BLAST "blastx -db /cvmfs/soft.mugqic/CentOS6/genomes/blast_db/uniprot_sprot.trinotate_v2.0.pep -max_target_seqs 1 -outfmt '6 std stitle'" 
blastx_trinity_uniprot.uniprot_sprot.trinotate_v2.0.pep.chunk_0000001.aa0c502502949b69d484967362fdb1a8.mugqic.done
chmod 755 $COMMAND
blastx_trinity_uniprot_2_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&     $COMMAND 
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE


if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem-per-cpu 4700M -c 40 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$blastx_trinity_uniprot_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$blastx_trinity_uniprot_2_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: blastx_trinity_uniprot_3_JOB_ID: blastx_trinity_uniprot.uniprot_sprot.trinotate_v2.0.pep.chunk_0000002
#-------------------------------------------------------------------------------
JOB_NAME=blastx_trinity_uniprot.uniprot_sprot.trinotate_v2.0.pep.chunk_0000002
JOB_DEPENDENCIES=$exonerate_fastasplit_1_JOB_ID
JOB_DONE=job_output/blastx_trinity_uniprot/blastx_trinity_uniprot.uniprot_sprot.trinotate_v2.0.pep.chunk_0000002.900e6a0d0ab53b40ae94fa7a9d9fe7c0.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'blastx_trinity_uniprot.uniprot_sprot.trinotate_v2.0.pep.chunk_0000002.900e6a0d0ab53b40ae94fa7a9d9fe7c0.mugqic.done' > $COMMAND
module purge && \
module load mugqic/perl/5.22.1 mugqic/mugqic_tools/2.8.2 mugqic/blast/2.3.0+ && \
mkdir -p blast && \
ln -s -f ../trinity_out_dir/Trinity.fasta_chunks/Trinity.fa_chunk_0000002 blast/query_Trinity_uniprot_sprot.trinotate_v2.0.pep_chunk_0000002.tsv && \
parallelBlast.pl \
-file blast/query_Trinity_uniprot_sprot.trinotate_v2.0.pep_chunk_0000002.tsv \
--OUT blast/blastx_Trinity_uniprot_sprot.trinotate_v2.0.pep_chunk_0000002.tsv \
-n 40 \
--BLAST "blastx -db /cvmfs/soft.mugqic/CentOS6/genomes/blast_db/uniprot_sprot.trinotate_v2.0.pep -max_target_seqs 1 -outfmt '6 std stitle'" 
blastx_trinity_uniprot.uniprot_sprot.trinotate_v2.0.pep.chunk_0000002.900e6a0d0ab53b40ae94fa7a9d9fe7c0.mugqic.done
chmod 755 $COMMAND
blastx_trinity_uniprot_3_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&     $COMMAND 
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE


if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem-per-cpu 4700M -c 40 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$blastx_trinity_uniprot_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$blastx_trinity_uniprot_3_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: blastx_trinity_uniprot_4_JOB_ID: blastx_trinity_uniprot.uniprot_sprot.trinotate_v2.0.pep.chunk_0000003
#-------------------------------------------------------------------------------
JOB_NAME=blastx_trinity_uniprot.uniprot_sprot.trinotate_v2.0.pep.chunk_0000003
JOB_DEPENDENCIES=$exonerate_fastasplit_1_JOB_ID
JOB_DONE=job_output/blastx_trinity_uniprot/blastx_trinity_uniprot.uniprot_sprot.trinotate_v2.0.pep.chunk_0000003.bb26f6f7a5d0fdfd8747d66b220ca4c3.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'blastx_trinity_uniprot.uniprot_sprot.trinotate_v2.0.pep.chunk_0000003.bb26f6f7a5d0fdfd8747d66b220ca4c3.mugqic.done' > $COMMAND
module purge && \
module load mugqic/perl/5.22.1 mugqic/mugqic_tools/2.8.2 mugqic/blast/2.3.0+ && \
mkdir -p blast && \
ln -s -f ../trinity_out_dir/Trinity.fasta_chunks/Trinity.fa_chunk_0000003 blast/query_Trinity_uniprot_sprot.trinotate_v2.0.pep_chunk_0000003.tsv && \
parallelBlast.pl \
-file blast/query_Trinity_uniprot_sprot.trinotate_v2.0.pep_chunk_0000003.tsv \
--OUT blast/blastx_Trinity_uniprot_sprot.trinotate_v2.0.pep_chunk_0000003.tsv \
-n 40 \
--BLAST "blastx -db /cvmfs/soft.mugqic/CentOS6/genomes/blast_db/uniprot_sprot.trinotate_v2.0.pep -max_target_seqs 1 -outfmt '6 std stitle'" 
blastx_trinity_uniprot.uniprot_sprot.trinotate_v2.0.pep.chunk_0000003.bb26f6f7a5d0fdfd8747d66b220ca4c3.mugqic.done
chmod 755 $COMMAND
blastx_trinity_uniprot_4_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&     $COMMAND 
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE


if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem-per-cpu 4700M -c 40 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$blastx_trinity_uniprot_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$blastx_trinity_uniprot_4_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: blastx_trinity_uniprot_5_JOB_ID: blastx_trinity_uniprot.uniprot_sprot.trinotate_v2.0.pep.chunk_0000004
#-------------------------------------------------------------------------------
JOB_NAME=blastx_trinity_uniprot.uniprot_sprot.trinotate_v2.0.pep.chunk_0000004
JOB_DEPENDENCIES=$exonerate_fastasplit_1_JOB_ID
JOB_DONE=job_output/blastx_trinity_uniprot/blastx_trinity_uniprot.uniprot_sprot.trinotate_v2.0.pep.chunk_0000004.2ca573639107fdc07b2eb0d62d22f709.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'blastx_trinity_uniprot.uniprot_sprot.trinotate_v2.0.pep.chunk_0000004.2ca573639107fdc07b2eb0d62d22f709.mugqic.done' > $COMMAND
module purge && \
module load mugqic/perl/5.22.1 mugqic/mugqic_tools/2.8.2 mugqic/blast/2.3.0+ && \
mkdir -p blast && \
ln -s -f ../trinity_out_dir/Trinity.fasta_chunks/Trinity.fa_chunk_0000004 blast/query_Trinity_uniprot_sprot.trinotate_v2.0.pep_chunk_0000004.tsv && \
parallelBlast.pl \
-file blast/query_Trinity_uniprot_sprot.trinotate_v2.0.pep_chunk_0000004.tsv \
--OUT blast/blastx_Trinity_uniprot_sprot.trinotate_v2.0.pep_chunk_0000004.tsv \
-n 40 \
--BLAST "blastx -db /cvmfs/soft.mugqic/CentOS6/genomes/blast_db/uniprot_sprot.trinotate_v2.0.pep -max_target_seqs 1 -outfmt '6 std stitle'" 
blastx_trinity_uniprot.uniprot_sprot.trinotate_v2.0.pep.chunk_0000004.2ca573639107fdc07b2eb0d62d22f709.mugqic.done
chmod 755 $COMMAND
blastx_trinity_uniprot_5_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&     $COMMAND 
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE


if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem-per-cpu 4700M -c 40 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$blastx_trinity_uniprot_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$blastx_trinity_uniprot_5_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: blastx_trinity_uniprot_6_JOB_ID: blastx_trinity_uniprot.uniprot_sprot.trinotate_v2.0.pep.chunk_0000005
#-------------------------------------------------------------------------------
JOB_NAME=blastx_trinity_uniprot.uniprot_sprot.trinotate_v2.0.pep.chunk_0000005
JOB_DEPENDENCIES=$exonerate_fastasplit_1_JOB_ID
JOB_DONE=job_output/blastx_trinity_uniprot/blastx_trinity_uniprot.uniprot_sprot.trinotate_v2.0.pep.chunk_0000005.8284b86f0a110d903f2c5ac5e1c3b164.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'blastx_trinity_uniprot.uniprot_sprot.trinotate_v2.0.pep.chunk_0000005.8284b86f0a110d903f2c5ac5e1c3b164.mugqic.done' > $COMMAND
module purge && \
module load mugqic/perl/5.22.1 mugqic/mugqic_tools/2.8.2 mugqic/blast/2.3.0+ && \
mkdir -p blast && \
ln -s -f ../trinity_out_dir/Trinity.fasta_chunks/Trinity.fa_chunk_0000005 blast/query_Trinity_uniprot_sprot.trinotate_v2.0.pep_chunk_0000005.tsv && \
parallelBlast.pl \
-file blast/query_Trinity_uniprot_sprot.trinotate_v2.0.pep_chunk_0000005.tsv \
--OUT blast/blastx_Trinity_uniprot_sprot.trinotate_v2.0.pep_chunk_0000005.tsv \
-n 40 \
--BLAST "blastx -db /cvmfs/soft.mugqic/CentOS6/genomes/blast_db/uniprot_sprot.trinotate_v2.0.pep -max_target_seqs 1 -outfmt '6 std stitle'" 
blastx_trinity_uniprot.uniprot_sprot.trinotate_v2.0.pep.chunk_0000005.8284b86f0a110d903f2c5ac5e1c3b164.mugqic.done
chmod 755 $COMMAND
blastx_trinity_uniprot_6_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&     $COMMAND 
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE


if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem-per-cpu 4700M -c 40 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$blastx_trinity_uniprot_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$blastx_trinity_uniprot_6_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: blastx_trinity_uniprot_7_JOB_ID: blastx_trinity_uniprot.uniprot_sprot.trinotate_v2.0.pep.chunk_0000006
#-------------------------------------------------------------------------------
JOB_NAME=blastx_trinity_uniprot.uniprot_sprot.trinotate_v2.0.pep.chunk_0000006
JOB_DEPENDENCIES=$exonerate_fastasplit_1_JOB_ID
JOB_DONE=job_output/blastx_trinity_uniprot/blastx_trinity_uniprot.uniprot_sprot.trinotate_v2.0.pep.chunk_0000006.56afbd4d35001fcacc4e46139e5cd594.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'blastx_trinity_uniprot.uniprot_sprot.trinotate_v2.0.pep.chunk_0000006.56afbd4d35001fcacc4e46139e5cd594.mugqic.done' > $COMMAND
module purge && \
module load mugqic/perl/5.22.1 mugqic/mugqic_tools/2.8.2 mugqic/blast/2.3.0+ && \
mkdir -p blast && \
ln -s -f ../trinity_out_dir/Trinity.fasta_chunks/Trinity.fa_chunk_0000006 blast/query_Trinity_uniprot_sprot.trinotate_v2.0.pep_chunk_0000006.tsv && \
parallelBlast.pl \
-file blast/query_Trinity_uniprot_sprot.trinotate_v2.0.pep_chunk_0000006.tsv \
--OUT blast/blastx_Trinity_uniprot_sprot.trinotate_v2.0.pep_chunk_0000006.tsv \
-n 40 \
--BLAST "blastx -db /cvmfs/soft.mugqic/CentOS6/genomes/blast_db/uniprot_sprot.trinotate_v2.0.pep -max_target_seqs 1 -outfmt '6 std stitle'" 
blastx_trinity_uniprot.uniprot_sprot.trinotate_v2.0.pep.chunk_0000006.56afbd4d35001fcacc4e46139e5cd594.mugqic.done
chmod 755 $COMMAND
blastx_trinity_uniprot_7_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&     $COMMAND 
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE


if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem-per-cpu 4700M -c 40 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$blastx_trinity_uniprot_7_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$blastx_trinity_uniprot_7_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: blastx_trinity_uniprot_8_JOB_ID: blastx_trinity_uniprot.uniprot_sprot.trinotate_v2.0.pep.chunk_0000007
#-------------------------------------------------------------------------------
JOB_NAME=blastx_trinity_uniprot.uniprot_sprot.trinotate_v2.0.pep.chunk_0000007
JOB_DEPENDENCIES=$exonerate_fastasplit_1_JOB_ID
JOB_DONE=job_output/blastx_trinity_uniprot/blastx_trinity_uniprot.uniprot_sprot.trinotate_v2.0.pep.chunk_0000007.8c0dab44ad888f282625c93f8d5b89ff.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'blastx_trinity_uniprot.uniprot_sprot.trinotate_v2.0.pep.chunk_0000007.8c0dab44ad888f282625c93f8d5b89ff.mugqic.done' > $COMMAND
module purge && \
module load mugqic/perl/5.22.1 mugqic/mugqic_tools/2.8.2 mugqic/blast/2.3.0+ && \
mkdir -p blast && \
ln -s -f ../trinity_out_dir/Trinity.fasta_chunks/Trinity.fa_chunk_0000007 blast/query_Trinity_uniprot_sprot.trinotate_v2.0.pep_chunk_0000007.tsv && \
parallelBlast.pl \
-file blast/query_Trinity_uniprot_sprot.trinotate_v2.0.pep_chunk_0000007.tsv \
--OUT blast/blastx_Trinity_uniprot_sprot.trinotate_v2.0.pep_chunk_0000007.tsv \
-n 40 \
--BLAST "blastx -db /cvmfs/soft.mugqic/CentOS6/genomes/blast_db/uniprot_sprot.trinotate_v2.0.pep -max_target_seqs 1 -outfmt '6 std stitle'" 
blastx_trinity_uniprot.uniprot_sprot.trinotate_v2.0.pep.chunk_0000007.8c0dab44ad888f282625c93f8d5b89ff.mugqic.done
chmod 755 $COMMAND
blastx_trinity_uniprot_8_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&     $COMMAND 
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE


if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem-per-cpu 4700M -c 40 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$blastx_trinity_uniprot_8_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$blastx_trinity_uniprot_8_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: blastx_trinity_uniprot_9_JOB_ID: blastx_trinity_uniprot.uniprot_sprot.trinotate_v2.0.pep.chunk_0000008
#-------------------------------------------------------------------------------
JOB_NAME=blastx_trinity_uniprot.uniprot_sprot.trinotate_v2.0.pep.chunk_0000008
JOB_DEPENDENCIES=$exonerate_fastasplit_1_JOB_ID
JOB_DONE=job_output/blastx_trinity_uniprot/blastx_trinity_uniprot.uniprot_sprot.trinotate_v2.0.pep.chunk_0000008.b10e6c027ba6e4c657af1fa735535418.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'blastx_trinity_uniprot.uniprot_sprot.trinotate_v2.0.pep.chunk_0000008.b10e6c027ba6e4c657af1fa735535418.mugqic.done' > $COMMAND
module purge && \
module load mugqic/perl/5.22.1 mugqic/mugqic_tools/2.8.2 mugqic/blast/2.3.0+ && \
mkdir -p blast && \
ln -s -f ../trinity_out_dir/Trinity.fasta_chunks/Trinity.fa_chunk_0000008 blast/query_Trinity_uniprot_sprot.trinotate_v2.0.pep_chunk_0000008.tsv && \
parallelBlast.pl \
-file blast/query_Trinity_uniprot_sprot.trinotate_v2.0.pep_chunk_0000008.tsv \
--OUT blast/blastx_Trinity_uniprot_sprot.trinotate_v2.0.pep_chunk_0000008.tsv \
-n 40 \
--BLAST "blastx -db /cvmfs/soft.mugqic/CentOS6/genomes/blast_db/uniprot_sprot.trinotate_v2.0.pep -max_target_seqs 1 -outfmt '6 std stitle'" 
blastx_trinity_uniprot.uniprot_sprot.trinotate_v2.0.pep.chunk_0000008.b10e6c027ba6e4c657af1fa735535418.mugqic.done
chmod 755 $COMMAND
blastx_trinity_uniprot_9_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&     $COMMAND 
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE


if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem-per-cpu 4700M -c 40 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$blastx_trinity_uniprot_9_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$blastx_trinity_uniprot_9_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: blastx_trinity_uniprot_10_JOB_ID: blastx_trinity_uniprot.uniprot_sprot.trinotate_v2.0.pep.chunk_0000009
#-------------------------------------------------------------------------------
JOB_NAME=blastx_trinity_uniprot.uniprot_sprot.trinotate_v2.0.pep.chunk_0000009
JOB_DEPENDENCIES=$exonerate_fastasplit_1_JOB_ID
JOB_DONE=job_output/blastx_trinity_uniprot/blastx_trinity_uniprot.uniprot_sprot.trinotate_v2.0.pep.chunk_0000009.048dc8c834406b0a194f6ed54dd72ad5.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'blastx_trinity_uniprot.uniprot_sprot.trinotate_v2.0.pep.chunk_0000009.048dc8c834406b0a194f6ed54dd72ad5.mugqic.done' > $COMMAND
module purge && \
module load mugqic/perl/5.22.1 mugqic/mugqic_tools/2.8.2 mugqic/blast/2.3.0+ && \
mkdir -p blast && \
ln -s -f ../trinity_out_dir/Trinity.fasta_chunks/Trinity.fa_chunk_0000009 blast/query_Trinity_uniprot_sprot.trinotate_v2.0.pep_chunk_0000009.tsv && \
parallelBlast.pl \
-file blast/query_Trinity_uniprot_sprot.trinotate_v2.0.pep_chunk_0000009.tsv \
--OUT blast/blastx_Trinity_uniprot_sprot.trinotate_v2.0.pep_chunk_0000009.tsv \
-n 40 \
--BLAST "blastx -db /cvmfs/soft.mugqic/CentOS6/genomes/blast_db/uniprot_sprot.trinotate_v2.0.pep -max_target_seqs 1 -outfmt '6 std stitle'" 
blastx_trinity_uniprot.uniprot_sprot.trinotate_v2.0.pep.chunk_0000009.048dc8c834406b0a194f6ed54dd72ad5.mugqic.done
chmod 755 $COMMAND
blastx_trinity_uniprot_10_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&     $COMMAND 
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE


if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem-per-cpu 4700M -c 40 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$blastx_trinity_uniprot_10_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$blastx_trinity_uniprot_10_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: blastx_trinity_uniprot_11_JOB_ID: blastx_trinity_uniprot.uniprot_sprot.trinotate_v2.0.pep.chunk_0000010
#-------------------------------------------------------------------------------
JOB_NAME=blastx_trinity_uniprot.uniprot_sprot.trinotate_v2.0.pep.chunk_0000010
JOB_DEPENDENCIES=$exonerate_fastasplit_1_JOB_ID
JOB_DONE=job_output/blastx_trinity_uniprot/blastx_trinity_uniprot.uniprot_sprot.trinotate_v2.0.pep.chunk_0000010.bb4ffecb2b0812be7621360bb8d36dea.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'blastx_trinity_uniprot.uniprot_sprot.trinotate_v2.0.pep.chunk_0000010.bb4ffecb2b0812be7621360bb8d36dea.mugqic.done' > $COMMAND
module purge && \
module load mugqic/perl/5.22.1 mugqic/mugqic_tools/2.8.2 mugqic/blast/2.3.0+ && \
mkdir -p blast && \
ln -s -f ../trinity_out_dir/Trinity.fasta_chunks/Trinity.fa_chunk_0000010 blast/query_Trinity_uniprot_sprot.trinotate_v2.0.pep_chunk_0000010.tsv && \
parallelBlast.pl \
-file blast/query_Trinity_uniprot_sprot.trinotate_v2.0.pep_chunk_0000010.tsv \
--OUT blast/blastx_Trinity_uniprot_sprot.trinotate_v2.0.pep_chunk_0000010.tsv \
-n 40 \
--BLAST "blastx -db /cvmfs/soft.mugqic/CentOS6/genomes/blast_db/uniprot_sprot.trinotate_v2.0.pep -max_target_seqs 1 -outfmt '6 std stitle'" 
blastx_trinity_uniprot.uniprot_sprot.trinotate_v2.0.pep.chunk_0000010.bb4ffecb2b0812be7621360bb8d36dea.mugqic.done
chmod 755 $COMMAND
blastx_trinity_uniprot_11_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&     $COMMAND 
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE


if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem-per-cpu 4700M -c 40 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$blastx_trinity_uniprot_11_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$blastx_trinity_uniprot_11_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: blastx_trinity_uniprot_12_JOB_ID: blastx_trinity_uniprot.uniprot_sprot.trinotate_v2.0.pep.chunk_0000011
#-------------------------------------------------------------------------------
JOB_NAME=blastx_trinity_uniprot.uniprot_sprot.trinotate_v2.0.pep.chunk_0000011
JOB_DEPENDENCIES=$exonerate_fastasplit_1_JOB_ID
JOB_DONE=job_output/blastx_trinity_uniprot/blastx_trinity_uniprot.uniprot_sprot.trinotate_v2.0.pep.chunk_0000011.da736fbc5dd8a11bf7a844d0a3e6645a.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'blastx_trinity_uniprot.uniprot_sprot.trinotate_v2.0.pep.chunk_0000011.da736fbc5dd8a11bf7a844d0a3e6645a.mugqic.done' > $COMMAND
module purge && \
module load mugqic/perl/5.22.1 mugqic/mugqic_tools/2.8.2 mugqic/blast/2.3.0+ && \
mkdir -p blast && \
ln -s -f ../trinity_out_dir/Trinity.fasta_chunks/Trinity.fa_chunk_0000011 blast/query_Trinity_uniprot_sprot.trinotate_v2.0.pep_chunk_0000011.tsv && \
parallelBlast.pl \
-file blast/query_Trinity_uniprot_sprot.trinotate_v2.0.pep_chunk_0000011.tsv \
--OUT blast/blastx_Trinity_uniprot_sprot.trinotate_v2.0.pep_chunk_0000011.tsv \
-n 40 \
--BLAST "blastx -db /cvmfs/soft.mugqic/CentOS6/genomes/blast_db/uniprot_sprot.trinotate_v2.0.pep -max_target_seqs 1 -outfmt '6 std stitle'" 
blastx_trinity_uniprot.uniprot_sprot.trinotate_v2.0.pep.chunk_0000011.da736fbc5dd8a11bf7a844d0a3e6645a.mugqic.done
chmod 755 $COMMAND
blastx_trinity_uniprot_12_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&     $COMMAND 
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE


if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem-per-cpu 4700M -c 40 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$blastx_trinity_uniprot_12_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$blastx_trinity_uniprot_12_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: blastx_trinity_uniprot_13_JOB_ID: blastx_trinity_uniprot.uniprot_sprot.trinotate_v2.0.pep.chunk_0000012
#-------------------------------------------------------------------------------
JOB_NAME=blastx_trinity_uniprot.uniprot_sprot.trinotate_v2.0.pep.chunk_0000012
JOB_DEPENDENCIES=$exonerate_fastasplit_1_JOB_ID
JOB_DONE=job_output/blastx_trinity_uniprot/blastx_trinity_uniprot.uniprot_sprot.trinotate_v2.0.pep.chunk_0000012.1d2b8c02a85e0e88aff8d1ec832dea74.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'blastx_trinity_uniprot.uniprot_sprot.trinotate_v2.0.pep.chunk_0000012.1d2b8c02a85e0e88aff8d1ec832dea74.mugqic.done' > $COMMAND
module purge && \
module load mugqic/perl/5.22.1 mugqic/mugqic_tools/2.8.2 mugqic/blast/2.3.0+ && \
mkdir -p blast && \
ln -s -f ../trinity_out_dir/Trinity.fasta_chunks/Trinity.fa_chunk_0000012 blast/query_Trinity_uniprot_sprot.trinotate_v2.0.pep_chunk_0000012.tsv && \
parallelBlast.pl \
-file blast/query_Trinity_uniprot_sprot.trinotate_v2.0.pep_chunk_0000012.tsv \
--OUT blast/blastx_Trinity_uniprot_sprot.trinotate_v2.0.pep_chunk_0000012.tsv \
-n 40 \
--BLAST "blastx -db /cvmfs/soft.mugqic/CentOS6/genomes/blast_db/uniprot_sprot.trinotate_v2.0.pep -max_target_seqs 1 -outfmt '6 std stitle'" 
blastx_trinity_uniprot.uniprot_sprot.trinotate_v2.0.pep.chunk_0000012.1d2b8c02a85e0e88aff8d1ec832dea74.mugqic.done
chmod 755 $COMMAND
blastx_trinity_uniprot_13_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&     $COMMAND 
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE


if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem-per-cpu 4700M -c 40 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$blastx_trinity_uniprot_13_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$blastx_trinity_uniprot_13_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: blastx_trinity_uniprot_14_JOB_ID: blastx_trinity_uniprot.uniprot_sprot.trinotate_v2.0.pep.chunk_0000013
#-------------------------------------------------------------------------------
JOB_NAME=blastx_trinity_uniprot.uniprot_sprot.trinotate_v2.0.pep.chunk_0000013
JOB_DEPENDENCIES=$exonerate_fastasplit_1_JOB_ID
JOB_DONE=job_output/blastx_trinity_uniprot/blastx_trinity_uniprot.uniprot_sprot.trinotate_v2.0.pep.chunk_0000013.e0724db073d2777547766fffa87f87bf.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'blastx_trinity_uniprot.uniprot_sprot.trinotate_v2.0.pep.chunk_0000013.e0724db073d2777547766fffa87f87bf.mugqic.done' > $COMMAND
module purge && \
module load mugqic/perl/5.22.1 mugqic/mugqic_tools/2.8.2 mugqic/blast/2.3.0+ && \
mkdir -p blast && \
ln -s -f ../trinity_out_dir/Trinity.fasta_chunks/Trinity.fa_chunk_0000013 blast/query_Trinity_uniprot_sprot.trinotate_v2.0.pep_chunk_0000013.tsv && \
parallelBlast.pl \
-file blast/query_Trinity_uniprot_sprot.trinotate_v2.0.pep_chunk_0000013.tsv \
--OUT blast/blastx_Trinity_uniprot_sprot.trinotate_v2.0.pep_chunk_0000013.tsv \
-n 40 \
--BLAST "blastx -db /cvmfs/soft.mugqic/CentOS6/genomes/blast_db/uniprot_sprot.trinotate_v2.0.pep -max_target_seqs 1 -outfmt '6 std stitle'" 
blastx_trinity_uniprot.uniprot_sprot.trinotate_v2.0.pep.chunk_0000013.e0724db073d2777547766fffa87f87bf.mugqic.done
chmod 755 $COMMAND
blastx_trinity_uniprot_14_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&     $COMMAND 
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE


if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem-per-cpu 4700M -c 40 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$blastx_trinity_uniprot_14_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$blastx_trinity_uniprot_14_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: blastx_trinity_uniprot_15_JOB_ID: blastx_trinity_uniprot.uniprot_sprot.trinotate_v2.0.pep.chunk_0000014
#-------------------------------------------------------------------------------
JOB_NAME=blastx_trinity_uniprot.uniprot_sprot.trinotate_v2.0.pep.chunk_0000014
JOB_DEPENDENCIES=$exonerate_fastasplit_1_JOB_ID
JOB_DONE=job_output/blastx_trinity_uniprot/blastx_trinity_uniprot.uniprot_sprot.trinotate_v2.0.pep.chunk_0000014.fba431c62b8369f2cd6de4e688a7487d.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'blastx_trinity_uniprot.uniprot_sprot.trinotate_v2.0.pep.chunk_0000014.fba431c62b8369f2cd6de4e688a7487d.mugqic.done' > $COMMAND
module purge && \
module load mugqic/perl/5.22.1 mugqic/mugqic_tools/2.8.2 mugqic/blast/2.3.0+ && \
mkdir -p blast && \
ln -s -f ../trinity_out_dir/Trinity.fasta_chunks/Trinity.fa_chunk_0000014 blast/query_Trinity_uniprot_sprot.trinotate_v2.0.pep_chunk_0000014.tsv && \
parallelBlast.pl \
-file blast/query_Trinity_uniprot_sprot.trinotate_v2.0.pep_chunk_0000014.tsv \
--OUT blast/blastx_Trinity_uniprot_sprot.trinotate_v2.0.pep_chunk_0000014.tsv \
-n 40 \
--BLAST "blastx -db /cvmfs/soft.mugqic/CentOS6/genomes/blast_db/uniprot_sprot.trinotate_v2.0.pep -max_target_seqs 1 -outfmt '6 std stitle'" 
blastx_trinity_uniprot.uniprot_sprot.trinotate_v2.0.pep.chunk_0000014.fba431c62b8369f2cd6de4e688a7487d.mugqic.done
chmod 755 $COMMAND
blastx_trinity_uniprot_15_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&     $COMMAND 
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE


if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem-per-cpu 4700M -c 40 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$blastx_trinity_uniprot_15_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$blastx_trinity_uniprot_15_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: blastx_trinity_uniprot_16_JOB_ID: blastx_trinity_uniprot.uniprot_sprot.trinotate_v2.0.pep.chunk_0000015
#-------------------------------------------------------------------------------
JOB_NAME=blastx_trinity_uniprot.uniprot_sprot.trinotate_v2.0.pep.chunk_0000015
JOB_DEPENDENCIES=$exonerate_fastasplit_1_JOB_ID
JOB_DONE=job_output/blastx_trinity_uniprot/blastx_trinity_uniprot.uniprot_sprot.trinotate_v2.0.pep.chunk_0000015.86bbe55f9604d757a635aedd84fdc1ea.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'blastx_trinity_uniprot.uniprot_sprot.trinotate_v2.0.pep.chunk_0000015.86bbe55f9604d757a635aedd84fdc1ea.mugqic.done' > $COMMAND
module purge && \
module load mugqic/perl/5.22.1 mugqic/mugqic_tools/2.8.2 mugqic/blast/2.3.0+ && \
mkdir -p blast && \
ln -s -f ../trinity_out_dir/Trinity.fasta_chunks/Trinity.fa_chunk_0000015 blast/query_Trinity_uniprot_sprot.trinotate_v2.0.pep_chunk_0000015.tsv && \
parallelBlast.pl \
-file blast/query_Trinity_uniprot_sprot.trinotate_v2.0.pep_chunk_0000015.tsv \
--OUT blast/blastx_Trinity_uniprot_sprot.trinotate_v2.0.pep_chunk_0000015.tsv \
-n 40 \
--BLAST "blastx -db /cvmfs/soft.mugqic/CentOS6/genomes/blast_db/uniprot_sprot.trinotate_v2.0.pep -max_target_seqs 1 -outfmt '6 std stitle'" 
blastx_trinity_uniprot.uniprot_sprot.trinotate_v2.0.pep.chunk_0000015.86bbe55f9604d757a635aedd84fdc1ea.mugqic.done
chmod 755 $COMMAND
blastx_trinity_uniprot_16_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&     $COMMAND 
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE


if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem-per-cpu 4700M -c 40 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$blastx_trinity_uniprot_16_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$blastx_trinity_uniprot_16_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: blastx_trinity_uniprot_17_JOB_ID: blastx_trinity_uniprot.uniprot_sprot.trinotate_v2.0.pep.chunk_0000016
#-------------------------------------------------------------------------------
JOB_NAME=blastx_trinity_uniprot.uniprot_sprot.trinotate_v2.0.pep.chunk_0000016
JOB_DEPENDENCIES=$exonerate_fastasplit_1_JOB_ID
JOB_DONE=job_output/blastx_trinity_uniprot/blastx_trinity_uniprot.uniprot_sprot.trinotate_v2.0.pep.chunk_0000016.634010a43bd00d4b02e60658f35866ff.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'blastx_trinity_uniprot.uniprot_sprot.trinotate_v2.0.pep.chunk_0000016.634010a43bd00d4b02e60658f35866ff.mugqic.done' > $COMMAND
module purge && \
module load mugqic/perl/5.22.1 mugqic/mugqic_tools/2.8.2 mugqic/blast/2.3.0+ && \
mkdir -p blast && \
ln -s -f ../trinity_out_dir/Trinity.fasta_chunks/Trinity.fa_chunk_0000016 blast/query_Trinity_uniprot_sprot.trinotate_v2.0.pep_chunk_0000016.tsv && \
parallelBlast.pl \
-file blast/query_Trinity_uniprot_sprot.trinotate_v2.0.pep_chunk_0000016.tsv \
--OUT blast/blastx_Trinity_uniprot_sprot.trinotate_v2.0.pep_chunk_0000016.tsv \
-n 40 \
--BLAST "blastx -db /cvmfs/soft.mugqic/CentOS6/genomes/blast_db/uniprot_sprot.trinotate_v2.0.pep -max_target_seqs 1 -outfmt '6 std stitle'" 
blastx_trinity_uniprot.uniprot_sprot.trinotate_v2.0.pep.chunk_0000016.634010a43bd00d4b02e60658f35866ff.mugqic.done
chmod 755 $COMMAND
blastx_trinity_uniprot_17_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&     $COMMAND 
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE


if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem-per-cpu 4700M -c 40 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$blastx_trinity_uniprot_17_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$blastx_trinity_uniprot_17_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: blastx_trinity_uniprot_18_JOB_ID: blastx_trinity_uniprot.uniprot_sprot.trinotate_v2.0.pep.chunk_0000017
#-------------------------------------------------------------------------------
JOB_NAME=blastx_trinity_uniprot.uniprot_sprot.trinotate_v2.0.pep.chunk_0000017
JOB_DEPENDENCIES=$exonerate_fastasplit_1_JOB_ID
JOB_DONE=job_output/blastx_trinity_uniprot/blastx_trinity_uniprot.uniprot_sprot.trinotate_v2.0.pep.chunk_0000017.2bc4917534b53d9cde59c56b5dde044a.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'blastx_trinity_uniprot.uniprot_sprot.trinotate_v2.0.pep.chunk_0000017.2bc4917534b53d9cde59c56b5dde044a.mugqic.done' > $COMMAND
module purge && \
module load mugqic/perl/5.22.1 mugqic/mugqic_tools/2.8.2 mugqic/blast/2.3.0+ && \
mkdir -p blast && \
ln -s -f ../trinity_out_dir/Trinity.fasta_chunks/Trinity.fa_chunk_0000017 blast/query_Trinity_uniprot_sprot.trinotate_v2.0.pep_chunk_0000017.tsv && \
parallelBlast.pl \
-file blast/query_Trinity_uniprot_sprot.trinotate_v2.0.pep_chunk_0000017.tsv \
--OUT blast/blastx_Trinity_uniprot_sprot.trinotate_v2.0.pep_chunk_0000017.tsv \
-n 40 \
--BLAST "blastx -db /cvmfs/soft.mugqic/CentOS6/genomes/blast_db/uniprot_sprot.trinotate_v2.0.pep -max_target_seqs 1 -outfmt '6 std stitle'" 
blastx_trinity_uniprot.uniprot_sprot.trinotate_v2.0.pep.chunk_0000017.2bc4917534b53d9cde59c56b5dde044a.mugqic.done
chmod 755 $COMMAND
blastx_trinity_uniprot_18_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&     $COMMAND 
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE


if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem-per-cpu 4700M -c 40 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$blastx_trinity_uniprot_18_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$blastx_trinity_uniprot_18_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: blastx_trinity_uniprot_19_JOB_ID: blastx_trinity_uniprot.uniprot_sprot.trinotate_v2.0.pep.chunk_0000018
#-------------------------------------------------------------------------------
JOB_NAME=blastx_trinity_uniprot.uniprot_sprot.trinotate_v2.0.pep.chunk_0000018
JOB_DEPENDENCIES=$exonerate_fastasplit_1_JOB_ID
JOB_DONE=job_output/blastx_trinity_uniprot/blastx_trinity_uniprot.uniprot_sprot.trinotate_v2.0.pep.chunk_0000018.fe923dfd029e3aa18b0b7ebe6ffe7814.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'blastx_trinity_uniprot.uniprot_sprot.trinotate_v2.0.pep.chunk_0000018.fe923dfd029e3aa18b0b7ebe6ffe7814.mugqic.done' > $COMMAND
module purge && \
module load mugqic/perl/5.22.1 mugqic/mugqic_tools/2.8.2 mugqic/blast/2.3.0+ && \
mkdir -p blast && \
ln -s -f ../trinity_out_dir/Trinity.fasta_chunks/Trinity.fa_chunk_0000018 blast/query_Trinity_uniprot_sprot.trinotate_v2.0.pep_chunk_0000018.tsv && \
parallelBlast.pl \
-file blast/query_Trinity_uniprot_sprot.trinotate_v2.0.pep_chunk_0000018.tsv \
--OUT blast/blastx_Trinity_uniprot_sprot.trinotate_v2.0.pep_chunk_0000018.tsv \
-n 40 \
--BLAST "blastx -db /cvmfs/soft.mugqic/CentOS6/genomes/blast_db/uniprot_sprot.trinotate_v2.0.pep -max_target_seqs 1 -outfmt '6 std stitle'" 
blastx_trinity_uniprot.uniprot_sprot.trinotate_v2.0.pep.chunk_0000018.fe923dfd029e3aa18b0b7ebe6ffe7814.mugqic.done
chmod 755 $COMMAND
blastx_trinity_uniprot_19_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&     $COMMAND 
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE


if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem-per-cpu 4700M -c 40 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$blastx_trinity_uniprot_19_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$blastx_trinity_uniprot_19_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: blastx_trinity_uniprot_20_JOB_ID: blastx_trinity_uniprot.uniprot_sprot.trinotate_v2.0.pep.chunk_0000019
#-------------------------------------------------------------------------------
JOB_NAME=blastx_trinity_uniprot.uniprot_sprot.trinotate_v2.0.pep.chunk_0000019
JOB_DEPENDENCIES=$exonerate_fastasplit_1_JOB_ID
JOB_DONE=job_output/blastx_trinity_uniprot/blastx_trinity_uniprot.uniprot_sprot.trinotate_v2.0.pep.chunk_0000019.d077b9cb3d068489a2f68018db144c3a.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'blastx_trinity_uniprot.uniprot_sprot.trinotate_v2.0.pep.chunk_0000019.d077b9cb3d068489a2f68018db144c3a.mugqic.done' > $COMMAND
module purge && \
module load mugqic/perl/5.22.1 mugqic/mugqic_tools/2.8.2 mugqic/blast/2.3.0+ && \
mkdir -p blast && \
ln -s -f ../trinity_out_dir/Trinity.fasta_chunks/Trinity.fa_chunk_0000019 blast/query_Trinity_uniprot_sprot.trinotate_v2.0.pep_chunk_0000019.tsv && \
parallelBlast.pl \
-file blast/query_Trinity_uniprot_sprot.trinotate_v2.0.pep_chunk_0000019.tsv \
--OUT blast/blastx_Trinity_uniprot_sprot.trinotate_v2.0.pep_chunk_0000019.tsv \
-n 40 \
--BLAST "blastx -db /cvmfs/soft.mugqic/CentOS6/genomes/blast_db/uniprot_sprot.trinotate_v2.0.pep -max_target_seqs 1 -outfmt '6 std stitle'" 
blastx_trinity_uniprot.uniprot_sprot.trinotate_v2.0.pep.chunk_0000019.d077b9cb3d068489a2f68018db144c3a.mugqic.done
chmod 755 $COMMAND
blastx_trinity_uniprot_20_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&     $COMMAND 
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE


if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem-per-cpu 4700M -c 40 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$blastx_trinity_uniprot_20_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$blastx_trinity_uniprot_20_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# STEP: blastx_trinity_uniprot_merge
#-------------------------------------------------------------------------------
STEP=blastx_trinity_uniprot_merge
mkdir -p $JOB_OUTPUT_DIR/$STEP

#-------------------------------------------------------------------------------
# JOB: blastx_trinity_uniprot_merge_1_JOB_ID: blastx_trinity_uniprot_sprot.trinotate_v2.0.pep_merge
#-------------------------------------------------------------------------------
JOB_NAME=blastx_trinity_uniprot_sprot.trinotate_v2.0.pep_merge
JOB_DEPENDENCIES=$blastx_trinity_uniprot_1_JOB_ID:$blastx_trinity_uniprot_2_JOB_ID:$blastx_trinity_uniprot_3_JOB_ID:$blastx_trinity_uniprot_4_JOB_ID:$blastx_trinity_uniprot_5_JOB_ID:$blastx_trinity_uniprot_6_JOB_ID:$blastx_trinity_uniprot_7_JOB_ID:$blastx_trinity_uniprot_8_JOB_ID:$blastx_trinity_uniprot_9_JOB_ID:$blastx_trinity_uniprot_10_JOB_ID:$blastx_trinity_uniprot_11_JOB_ID:$blastx_trinity_uniprot_12_JOB_ID:$blastx_trinity_uniprot_13_JOB_ID:$blastx_trinity_uniprot_14_JOB_ID:$blastx_trinity_uniprot_15_JOB_ID:$blastx_trinity_uniprot_16_JOB_ID:$blastx_trinity_uniprot_17_JOB_ID:$blastx_trinity_uniprot_18_JOB_ID:$blastx_trinity_uniprot_19_JOB_ID:$blastx_trinity_uniprot_20_JOB_ID
JOB_DONE=job_output/blastx_trinity_uniprot_merge/blastx_trinity_uniprot_sprot.trinotate_v2.0.pep_merge.e081d08ef7d5216625f3ac70d63204e2.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'blastx_trinity_uniprot_sprot.trinotate_v2.0.pep_merge.e081d08ef7d5216625f3ac70d63204e2.mugqic.done' > $COMMAND
cat \
  blast/blastx_Trinity_uniprot_sprot.trinotate_v2.0.pep_chunk_0000000.tsv \
  blast/blastx_Trinity_uniprot_sprot.trinotate_v2.0.pep_chunk_0000001.tsv \
  blast/blastx_Trinity_uniprot_sprot.trinotate_v2.0.pep_chunk_0000002.tsv \
  blast/blastx_Trinity_uniprot_sprot.trinotate_v2.0.pep_chunk_0000003.tsv \
  blast/blastx_Trinity_uniprot_sprot.trinotate_v2.0.pep_chunk_0000004.tsv \
  blast/blastx_Trinity_uniprot_sprot.trinotate_v2.0.pep_chunk_0000005.tsv \
  blast/blastx_Trinity_uniprot_sprot.trinotate_v2.0.pep_chunk_0000006.tsv \
  blast/blastx_Trinity_uniprot_sprot.trinotate_v2.0.pep_chunk_0000007.tsv \
  blast/blastx_Trinity_uniprot_sprot.trinotate_v2.0.pep_chunk_0000008.tsv \
  blast/blastx_Trinity_uniprot_sprot.trinotate_v2.0.pep_chunk_0000009.tsv \
  blast/blastx_Trinity_uniprot_sprot.trinotate_v2.0.pep_chunk_0000010.tsv \
  blast/blastx_Trinity_uniprot_sprot.trinotate_v2.0.pep_chunk_0000011.tsv \
  blast/blastx_Trinity_uniprot_sprot.trinotate_v2.0.pep_chunk_0000012.tsv \
  blast/blastx_Trinity_uniprot_sprot.trinotate_v2.0.pep_chunk_0000013.tsv \
  blast/blastx_Trinity_uniprot_sprot.trinotate_v2.0.pep_chunk_0000014.tsv \
  blast/blastx_Trinity_uniprot_sprot.trinotate_v2.0.pep_chunk_0000015.tsv \
  blast/blastx_Trinity_uniprot_sprot.trinotate_v2.0.pep_chunk_0000016.tsv \
  blast/blastx_Trinity_uniprot_sprot.trinotate_v2.0.pep_chunk_0000017.tsv \
  blast/blastx_Trinity_uniprot_sprot.trinotate_v2.0.pep_chunk_0000018.tsv \
  blast/blastx_Trinity_uniprot_sprot.trinotate_v2.0.pep_chunk_0000019.tsv \
  > blast/blastx_Trinity_uniprot_sprot.trinotate_v2.0.pep.tsv && \
zip -j blast/blastx_Trinity_uniprot_sprot.trinotate_v2.0.pep.tsv.zip blast/blastx_Trinity_uniprot_sprot.trinotate_v2.0.pep.tsv
blastx_trinity_uniprot_sprot.trinotate_v2.0.pep_merge.e081d08ef7d5216625f3ac70d63204e2.mugqic.done
chmod 755 $COMMAND
blastx_trinity_uniprot_merge_1_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&     $COMMAND 
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE


if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem-per-cpu 4700M -c 1 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$blastx_trinity_uniprot_merge_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$blastx_trinity_uniprot_merge_1_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: blastx_trinity_uniprot_merge_2_JOB_ID: blastx_trinity_uniprot_merge_report
#-------------------------------------------------------------------------------
JOB_NAME=blastx_trinity_uniprot_merge_report
JOB_DEPENDENCIES=$blastx_trinity_uniprot_merge_1_JOB_ID
JOB_DONE=job_output/blastx_trinity_uniprot_merge/blastx_trinity_uniprot_merge_report.908bcced330ffaca21522bef995e6bd0.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'blastx_trinity_uniprot_merge_report.908bcced330ffaca21522bef995e6bd0.mugqic.done' > $COMMAND
module purge && \
module load mugqic/pandoc/2.16.1 && \
mkdir -p report && \
cp blast/blastx_Trinity_uniprot_sprot.trinotate_v2.0.pep.tsv.zip  report/ && \
pandoc --to=markdown \
  --template /cvmfs/soft.mugqic/root/software/genpipes/genpipes-4.3.1/bfx/report/RnaSeqDeNovoAssembly.blastx_trinity_uniprot_merge.md \
  --variable blast_db="uniprot_sprot.trinotate_v2.0.pep" \
  /cvmfs/soft.mugqic/root/software/genpipes/genpipes-4.3.1/bfx/report/RnaSeqDeNovoAssembly.blastx_trinity_uniprot_merge.md \
  > report/RnaSeqDeNovoAssembly.blastx_trinity_uniprot_merge.md
blastx_trinity_uniprot_merge_report.908bcced330ffaca21522bef995e6bd0.mugqic.done
chmod 755 $COMMAND
blastx_trinity_uniprot_merge_2_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&     $COMMAND 
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE


if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem-per-cpu 4700M -c 1 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$blastx_trinity_uniprot_merge_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$blastx_trinity_uniprot_merge_2_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# STEP: transdecoder
#-------------------------------------------------------------------------------
STEP=transdecoder
mkdir -p $JOB_OUTPUT_DIR/$STEP

#-------------------------------------------------------------------------------
# JOB: transdecoder_1_JOB_ID: transdecoder
#-------------------------------------------------------------------------------
JOB_NAME=transdecoder
JOB_DEPENDENCIES=$trinity_1_JOB_ID
JOB_DONE=job_output/transdecoder/transdecoder.cbfc06dda51fe755125f69addb1bbe62.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'transdecoder.cbfc06dda51fe755125f69addb1bbe62.mugqic.done' > $COMMAND
module purge && \
module load mugqic/perl/5.22.1 mugqic/blast/2.3.0+ mugqic/hmmer/3.1b2 mugqic/trinotate/2.0.2 mugqic/TransDecoder/2.0.1 && \
mkdir -p trinotate/transdecoder && \
cd trinotate/transdecoder && \
TransDecoder.LongOrfs -S \
-t ../../trinity_out_dir/Trinity.fasta && \
blastp \
-query Trinity.fasta.transdecoder_dir/longest_orfs.pep \
-db /cvmfs/soft.mugqic/CentOS6/genomes/blast_db/uniprot_sprot.trinotate_v2.0.pep \
-max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads 40 \
> Trinity.fasta.transdecoder_dir/longest_orfs.blastp.outfmt6 && \
hmmscan --cpu 40 \
--domtblout Trinity.fasta.transdecoder_dir/longest_orfs.pfam.domtblout \
/cvmfs/soft.mugqic/CentOS6/genomes/pfam_db/Pfam-A.hmm \
Trinity.fasta.transdecoder_dir/longest_orfs.pep && \
TransDecoder.Predict \
-t ../../trinity_out_dir/Trinity.fasta \
--retain_pfam_hits Trinity.fasta.transdecoder_dir/longest_orfs.pfam.domtblout \
--retain_blastp_hits Trinity.fasta.transdecoder_dir/longest_orfs.blastp.outfmt6 && \
cd ../..
transdecoder.cbfc06dda51fe755125f69addb1bbe62.mugqic.done
chmod 755 $COMMAND
transdecoder_1_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&     $COMMAND 
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE


if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem-per-cpu 4700M -c 40 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$transdecoder_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$transdecoder_1_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# STEP: hmmer
#-------------------------------------------------------------------------------
STEP=hmmer
mkdir -p $JOB_OUTPUT_DIR/$STEP

#-------------------------------------------------------------------------------
# JOB: hmmer_1_JOB_ID: hmmer
#-------------------------------------------------------------------------------
JOB_NAME=hmmer
JOB_DEPENDENCIES=$transdecoder_1_JOB_ID
JOB_DONE=job_output/hmmer/hmmer.0e83c0521eebc79ff4335296b5411b13.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'hmmer.0e83c0521eebc79ff4335296b5411b13.mugqic.done' > $COMMAND
module purge && \
module load mugqic/hmmer/3.1b2 && \
hmmscan --cpu 40 \
--domtblout trinotate/transdecoder/Trinity.fasta.transdecoder.pfam \
/cvmfs/soft.mugqic/CentOS6/genomes/pfam_db/Pfam-A.hmm \
trinotate/transdecoder/Trinity.fasta.transdecoder.pep
hmmer.0e83c0521eebc79ff4335296b5411b13.mugqic.done
chmod 755 $COMMAND
hmmer_1_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&     $COMMAND 
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE


if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem-per-cpu 4700M -c 40 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$hmmer_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$hmmer_1_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# STEP: rnammer_transcriptome
#-------------------------------------------------------------------------------
STEP=rnammer_transcriptome
mkdir -p $JOB_OUTPUT_DIR/$STEP

#-------------------------------------------------------------------------------
# JOB: rnammer_transcriptome_1_JOB_ID: rnammer_transcriptome
#-------------------------------------------------------------------------------
JOB_NAME=rnammer_transcriptome
JOB_DEPENDENCIES=$trinity_1_JOB_ID
JOB_DONE=job_output/rnammer_transcriptome/rnammer_transcriptome.58b905516a0c1a3f2beeeb3a50673c48.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'rnammer_transcriptome.58b905516a0c1a3f2beeeb3a50673c48.mugqic.done' > $COMMAND
module purge && \
module load mugqic/perl/5.22.1 mugqic/hmmer/2.3.2 mugqic/rnammer/1.2 mugqic/trinity/2.0.4_patch mugqic/trinotate/2.0.2 && \
mkdir -p trinotate/rnammer && \
cd trinotate/rnammer && \
$TRINOTATE_HOME/util/rnammer_support/RnammerTranscriptome.pl  \
--transcriptome ../../trinity_out_dir/Trinity.fasta \
--path_to_rnammer `which rnammer` && \
cd ../..
rnammer_transcriptome.58b905516a0c1a3f2beeeb3a50673c48.mugqic.done
chmod 755 $COMMAND
rnammer_transcriptome_1_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&     $COMMAND 
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE


if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem 20G -c 2 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$rnammer_transcriptome_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$rnammer_transcriptome_1_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# STEP: blastp_transdecoder_uniprot
#-------------------------------------------------------------------------------
STEP=blastp_transdecoder_uniprot
mkdir -p $JOB_OUTPUT_DIR/$STEP

#-------------------------------------------------------------------------------
# JOB: blastp_transdecoder_uniprot_1_JOB_ID: blastp_transdecoder_uniprot.uniprot_sprot.trinotate_v2.0.pep
#-------------------------------------------------------------------------------
JOB_NAME=blastp_transdecoder_uniprot.uniprot_sprot.trinotate_v2.0.pep
JOB_DEPENDENCIES=$transdecoder_1_JOB_ID
JOB_DONE=job_output/blastp_transdecoder_uniprot/blastp_transdecoder_uniprot.uniprot_sprot.trinotate_v2.0.pep.769b5ccdfc22c8943a5f8ba9a6976268.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'blastp_transdecoder_uniprot.uniprot_sprot.trinotate_v2.0.pep.769b5ccdfc22c8943a5f8ba9a6976268.mugqic.done' > $COMMAND
module purge && \
module load mugqic/perl/5.22.1 mugqic/mugqic_tools/2.8.2 mugqic/blast/2.3.0+ && \
mkdir -p trinotate/blastp && \
ln -s -f ../transdecoder/Trinity.fasta.transdecoder.pep trinotate/blastp/Trinity.fasta.transdecoder.pep_uniprot_sprot.trinotate_v2.0.pep.tsv && \
parallelBlast.pl \
-file trinotate/blastp/Trinity.fasta.transdecoder.pep_uniprot_sprot.trinotate_v2.0.pep.tsv \
--OUT trinotate/blastp/blastp_Trinity.fasta.transdecoder.pep_uniprot_sprot.trinotate_v2.0.pep.tsv \
-n 40 \
--BLAST "blastp -db /cvmfs/soft.mugqic/CentOS6/genomes/blast_db/uniprot_sprot.trinotate_v2.0.pep -max_target_seqs 1 -outfmt '6 std stitle'" 
blastp_transdecoder_uniprot.uniprot_sprot.trinotate_v2.0.pep.769b5ccdfc22c8943a5f8ba9a6976268.mugqic.done
chmod 755 $COMMAND
blastp_transdecoder_uniprot_1_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&     $COMMAND 
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE


if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem-per-cpu 4700M -c 40 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$blastp_transdecoder_uniprot_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$blastp_transdecoder_uniprot_1_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# STEP: signalp
#-------------------------------------------------------------------------------
STEP=signalp
mkdir -p $JOB_OUTPUT_DIR/$STEP

#-------------------------------------------------------------------------------
# JOB: signalp_1_JOB_ID: signalp
#-------------------------------------------------------------------------------
JOB_NAME=signalp
JOB_DEPENDENCIES=$transdecoder_1_JOB_ID
JOB_DONE=job_output/signalp/signalp.0dc676527682310f762feca08a4c17b6.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'signalp.0dc676527682310f762feca08a4c17b6.mugqic.done' > $COMMAND
module purge && \
module load mugqic/perl/5.22.1 mugqic/signalp/4.1 && \
signalp -f short  \
-T trinotate/signalp \
-n trinotate/signalp/signalp.gff \
trinotate/transdecoder/Trinity.fasta.transdecoder.pep
signalp.0dc676527682310f762feca08a4c17b6.mugqic.done
chmod 755 $COMMAND
signalp_1_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&     $COMMAND 
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE


if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem 20G -c 2 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$signalp_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$signalp_1_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# STEP: tmhmm
#-------------------------------------------------------------------------------
STEP=tmhmm
mkdir -p $JOB_OUTPUT_DIR/$STEP

#-------------------------------------------------------------------------------
# JOB: tmhmm_1_JOB_ID: tmhmm
#-------------------------------------------------------------------------------
JOB_NAME=tmhmm
JOB_DEPENDENCIES=$transdecoder_1_JOB_ID
JOB_DONE=job_output/tmhmm/tmhmm.739c90aab4806c006afb43fa76b17cf2.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'tmhmm.739c90aab4806c006afb43fa76b17cf2.mugqic.done' > $COMMAND
module purge && \
module load mugqic/perl/5.22.1 mugqic/tmhmm/2.0c && \
mkdir -p trinotate/tmhmm && \
tmhmm --short \
< trinotate/transdecoder/Trinity.fasta.transdecoder.pep \
> trinotate/tmhmm/tmhmm.out
tmhmm.739c90aab4806c006afb43fa76b17cf2.mugqic.done
chmod 755 $COMMAND
tmhmm_1_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&     $COMMAND 
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE


if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem-per-cpu 4700M -c 1 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$tmhmm_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$tmhmm_1_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# STEP: trinotate
#-------------------------------------------------------------------------------
STEP=trinotate
mkdir -p $JOB_OUTPUT_DIR/$STEP

#-------------------------------------------------------------------------------
# JOB: trinotate_1_JOB_ID: trinotate
#-------------------------------------------------------------------------------
JOB_NAME=trinotate
JOB_DEPENDENCIES=$trinity_1_JOB_ID:$blastx_trinity_uniprot_merge_1_JOB_ID:$transdecoder_1_JOB_ID:$hmmer_1_JOB_ID:$rnammer_transcriptome_1_JOB_ID:$blastp_transdecoder_uniprot_1_JOB_ID:$signalp_1_JOB_ID:$tmhmm_1_JOB_ID
JOB_DONE=job_output/trinotate/trinotate.9240debe20ba2016c470bfb46d3f0fad.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'trinotate.9240debe20ba2016c470bfb46d3f0fad.mugqic.done' > $COMMAND
module purge && \
module load mugqic/perl/5.22.1 mugqic/trinity/2.0.4_patch mugqic/trinotate/2.0.2 && \
mkdir -p trinotate && \
cp $TRINOTATE_SQLITE trinotate/Trinotate.sqlite && \
$TRINITY_HOME/util/support_scripts/get_Trinity_gene_to_trans_map.pl \
trinity_out_dir/Trinity.fasta \
> trinity_out_dir/Trinity.fasta.gene_trans_map && \
Trinotate trinotate/Trinotate.sqlite init \
--gene_trans_map trinity_out_dir/Trinity.fasta.gene_trans_map \
--transcript trinity_out_dir/Trinity.fasta \
--transdecoder_pep trinotate/transdecoder/Trinity.fasta.transdecoder.pep && \
Trinotate trinotate/Trinotate.sqlite LOAD_swissprot_blastx blast/blastx_Trinity_uniprot_sprot.trinotate_v2.0.pep.tsv && \
Trinotate trinotate/Trinotate.sqlite LOAD_swissprot_blastp trinotate/blastp/blastp_Trinity.fasta.transdecoder.pep_uniprot_sprot.trinotate_v2.0.pep.tsv && \
Trinotate trinotate/Trinotate.sqlite LOAD_pfam trinotate/transdecoder/Trinity.fasta.transdecoder.pfam && \
Trinotate trinotate/Trinotate.sqlite LOAD_tmhmm trinotate/tmhmm/tmhmm.out && \
Trinotate trinotate/Trinotate.sqlite LOAD_signalp trinotate/signalp/signalp.gff && \
Trinotate trinotate/Trinotate.sqlite LOAD_rnammer trinotate/rnammer/Trinity.fasta.rnammer.gff && \
Trinotate trinotate/Trinotate.sqlite report -E 1e-5 --pfam_cutoff DNC \
> trinotate/trinotate_annotation_report.tsv
trinotate.9240debe20ba2016c470bfb46d3f0fad.mugqic.done
chmod 755 $COMMAND
trinotate_1_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&     $COMMAND 
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE


if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem-per-cpu 4700M -c 1 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$trinotate_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$trinotate_1_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: trinotate_2_JOB_ID: trinotate_report
#-------------------------------------------------------------------------------
JOB_NAME=trinotate_report
JOB_DEPENDENCIES=$trinotate_1_JOB_ID
JOB_DONE=job_output/trinotate/trinotate_report.ef03c9696ae7792cd07b7445adf8deca.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'trinotate_report.ef03c9696ae7792cd07b7445adf8deca.mugqic.done' > $COMMAND
module purge && \
module load mugqic/R_Bioconductor/3.6.0_3.9 mugqic/pandoc/2.16.1 && \
R --no-save --no-restore <<-'EOF'
report_dir="report"; source_dir="trinotate";
input_rmarkdown_file = '/cvmfs/soft.mugqic/root/software/genpipes/genpipes-4.3.1/bfx/report/RnaSeqDeNovoAssembly.trinotate.Rmd'
render_output_dir    = 'report'
rmarkdown_file       = basename(input_rmarkdown_file) # honoring a different WD that location of Rmd file in knitr is problematic
file.copy(from = input_rmarkdown_file, to = rmarkdown_file, overwrite = T)
rmarkdown::render(input = rmarkdown_file, output_format = c("html_document","md_document"), output_dir = render_output_dir  )
file.remove(rmarkdown_file)
EOF
trinotate_report.ef03c9696ae7792cd07b7445adf8deca.mugqic.done
chmod 755 $COMMAND
trinotate_2_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&     $COMMAND 
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE


if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem-per-cpu 4700M -c 1 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$trinotate_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$trinotate_2_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# STEP: align_and_estimate_abundance_prep_reference
#-------------------------------------------------------------------------------
STEP=align_and_estimate_abundance_prep_reference
mkdir -p $JOB_OUTPUT_DIR/$STEP

#-------------------------------------------------------------------------------
# JOB: align_and_estimate_abundance_prep_reference_1_JOB_ID: align_and_estimate_abundance_prep_reference
#-------------------------------------------------------------------------------
JOB_NAME=align_and_estimate_abundance_prep_reference
JOB_DEPENDENCIES=$trinity_1_JOB_ID
JOB_DONE=job_output/align_and_estimate_abundance_prep_reference/align_and_estimate_abundance_prep_reference.f0509b0227cfd0a289530d94fb1e6ce4.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'align_and_estimate_abundance_prep_reference.f0509b0227cfd0a289530d94fb1e6ce4.mugqic.done' > $COMMAND
module purge && \
module load mugqic/perl/5.22.1 mugqic/bowtie/1.0.0 mugqic/samtools/1.12 mugqic/trinity/2.0.4_patch && \
align_and_estimate_abundance.pl \
  --transcripts trinity_out_dir/Trinity.fasta \
  --seqType fa \
  --est_method RSEM \
  --aln_method bowtie \
  --trinity_mode \
  --output_dir trinity_out_dir \
  --prep_reference
align_and_estimate_abundance_prep_reference.f0509b0227cfd0a289530d94fb1e6ce4.mugqic.done
chmod 755 $COMMAND
align_and_estimate_abundance_prep_reference_1_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&     $COMMAND 
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE


if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem-per-cpu 4700M -c 1 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$align_and_estimate_abundance_prep_reference_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$align_and_estimate_abundance_prep_reference_1_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# STEP: align_and_estimate_abundance
#-------------------------------------------------------------------------------
STEP=align_and_estimate_abundance
mkdir -p $JOB_OUTPUT_DIR/$STEP

#-------------------------------------------------------------------------------
# JOB: align_and_estimate_abundance_1_JOB_ID: align_and_estimate_abundance.GM12878_Rep1
#-------------------------------------------------------------------------------
JOB_NAME=align_and_estimate_abundance.GM12878_Rep1
JOB_DEPENDENCIES=$trimmomatic_1_JOB_ID:$trinity_1_JOB_ID:$align_and_estimate_abundance_prep_reference_1_JOB_ID
JOB_DONE=job_output/align_and_estimate_abundance/align_and_estimate_abundance.GM12878_Rep1.ae633223b71083cdeee24653772ce5bd.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'align_and_estimate_abundance.GM12878_Rep1.ae633223b71083cdeee24653772ce5bd.mugqic.done' > $COMMAND
module purge && \
module load mugqic/perl/5.22.1 mugqic/bowtie/1.0.0 mugqic/samtools/1.12 mugqic/trinity/2.0.4_patch && \
align_and_estimate_abundance.pl --SS_lib_type RF \
  --transcripts trinity_out_dir/Trinity.fasta \
  --seqType fq \
  --est_method RSEM \
  --aln_method bowtie \
  --trinity_mode \
  --output_prefix GM12878_Rep1 \
  --output_dir align_and_estimate_abundance/GM12878_Rep1 \
  --thread_count 5 \
  --left trim/GM12878_Rep1/GM12878_Rep1.trim.pair1.fastq.gz \
  --right trim/GM12878_Rep1/GM12878_Rep1.trim.pair2.fastq.gz
align_and_estimate_abundance.GM12878_Rep1.ae633223b71083cdeee24653772ce5bd.mugqic.done
chmod 755 $COMMAND
align_and_estimate_abundance_1_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&     $COMMAND 
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE


if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem-per-cpu 4700M -c 5 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$align_and_estimate_abundance_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$align_and_estimate_abundance_1_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: align_and_estimate_abundance_2_JOB_ID: align_and_estimate_abundance.GM12878_Rep2
#-------------------------------------------------------------------------------
JOB_NAME=align_and_estimate_abundance.GM12878_Rep2
JOB_DEPENDENCIES=$trimmomatic_2_JOB_ID:$trinity_1_JOB_ID:$align_and_estimate_abundance_prep_reference_1_JOB_ID
JOB_DONE=job_output/align_and_estimate_abundance/align_and_estimate_abundance.GM12878_Rep2.d2e86015bf3828248dfceeeb699e48d6.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'align_and_estimate_abundance.GM12878_Rep2.d2e86015bf3828248dfceeeb699e48d6.mugqic.done' > $COMMAND
module purge && \
module load mugqic/perl/5.22.1 mugqic/bowtie/1.0.0 mugqic/samtools/1.12 mugqic/trinity/2.0.4_patch && \
align_and_estimate_abundance.pl --SS_lib_type RF \
  --transcripts trinity_out_dir/Trinity.fasta \
  --seqType fq \
  --est_method RSEM \
  --aln_method bowtie \
  --trinity_mode \
  --output_prefix GM12878_Rep2 \
  --output_dir align_and_estimate_abundance/GM12878_Rep2 \
  --thread_count 5 \
  --left trim/GM12878_Rep2/GM12878_Rep2.trim.pair1.fastq.gz \
  --right trim/GM12878_Rep2/GM12878_Rep2.trim.pair2.fastq.gz
align_and_estimate_abundance.GM12878_Rep2.d2e86015bf3828248dfceeeb699e48d6.mugqic.done
chmod 755 $COMMAND
align_and_estimate_abundance_2_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&     $COMMAND 
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE


if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem-per-cpu 4700M -c 5 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$align_and_estimate_abundance_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$align_and_estimate_abundance_2_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: align_and_estimate_abundance_3_JOB_ID: align_and_estimate_abundance.H1ESC_Rep1
#-------------------------------------------------------------------------------
JOB_NAME=align_and_estimate_abundance.H1ESC_Rep1
JOB_DEPENDENCIES=$trimmomatic_3_JOB_ID:$trinity_1_JOB_ID:$align_and_estimate_abundance_prep_reference_1_JOB_ID
JOB_DONE=job_output/align_and_estimate_abundance/align_and_estimate_abundance.H1ESC_Rep1.fa4eeaba8887eb654a81deab763116ba.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'align_and_estimate_abundance.H1ESC_Rep1.fa4eeaba8887eb654a81deab763116ba.mugqic.done' > $COMMAND
module purge && \
module load mugqic/perl/5.22.1 mugqic/bowtie/1.0.0 mugqic/samtools/1.12 mugqic/trinity/2.0.4_patch && \
align_and_estimate_abundance.pl --SS_lib_type RF \
  --transcripts trinity_out_dir/Trinity.fasta \
  --seqType fq \
  --est_method RSEM \
  --aln_method bowtie \
  --trinity_mode \
  --output_prefix H1ESC_Rep1 \
  --output_dir align_and_estimate_abundance/H1ESC_Rep1 \
  --thread_count 5 \
  --left trim/H1ESC_Rep1/H1ESC_Rep1.trim.pair1.fastq.gz \
  --right trim/H1ESC_Rep1/H1ESC_Rep1.trim.pair2.fastq.gz
align_and_estimate_abundance.H1ESC_Rep1.fa4eeaba8887eb654a81deab763116ba.mugqic.done
chmod 755 $COMMAND
align_and_estimate_abundance_3_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&     $COMMAND 
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE


if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem-per-cpu 4700M -c 5 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$align_and_estimate_abundance_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$align_and_estimate_abundance_3_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: align_and_estimate_abundance_4_JOB_ID: align_and_estimate_abundance.H1ESC_Rep2
#-------------------------------------------------------------------------------
JOB_NAME=align_and_estimate_abundance.H1ESC_Rep2
JOB_DEPENDENCIES=$trimmomatic_4_JOB_ID:$trinity_1_JOB_ID:$align_and_estimate_abundance_prep_reference_1_JOB_ID
JOB_DONE=job_output/align_and_estimate_abundance/align_and_estimate_abundance.H1ESC_Rep2.c95e1e7cf89064240d52f61507b46506.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'align_and_estimate_abundance.H1ESC_Rep2.c95e1e7cf89064240d52f61507b46506.mugqic.done' > $COMMAND
module purge && \
module load mugqic/perl/5.22.1 mugqic/bowtie/1.0.0 mugqic/samtools/1.12 mugqic/trinity/2.0.4_patch && \
align_and_estimate_abundance.pl --SS_lib_type RF \
  --transcripts trinity_out_dir/Trinity.fasta \
  --seqType fq \
  --est_method RSEM \
  --aln_method bowtie \
  --trinity_mode \
  --output_prefix H1ESC_Rep2 \
  --output_dir align_and_estimate_abundance/H1ESC_Rep2 \
  --thread_count 5 \
  --left trim/H1ESC_Rep2/H1ESC_Rep2.trim.pair1.fastq.gz \
  --right trim/H1ESC_Rep2/H1ESC_Rep2.trim.pair2.fastq.gz
align_and_estimate_abundance.H1ESC_Rep2.c95e1e7cf89064240d52f61507b46506.mugqic.done
chmod 755 $COMMAND
align_and_estimate_abundance_4_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&     $COMMAND 
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE


if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem-per-cpu 4700M -c 5 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$align_and_estimate_abundance_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$align_and_estimate_abundance_4_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: align_and_estimate_abundance_5_JOB_ID: align_and_estimate_abundance.genes
#-------------------------------------------------------------------------------
JOB_NAME=align_and_estimate_abundance.genes
JOB_DEPENDENCIES=$align_and_estimate_abundance_1_JOB_ID:$align_and_estimate_abundance_2_JOB_ID:$align_and_estimate_abundance_3_JOB_ID:$align_and_estimate_abundance_4_JOB_ID
JOB_DONE=job_output/align_and_estimate_abundance/align_and_estimate_abundance.genes.b3c4c1fd6af7e4eaf50b3947423ed48a.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'align_and_estimate_abundance.genes.b3c4c1fd6af7e4eaf50b3947423ed48a.mugqic.done' > $COMMAND
module purge && \
module load mugqic/perl/5.22.1 mugqic/trinity/2.0.4_patch mugqic/R_Bioconductor/3.6.0_3.9 && \
mkdir -p differential_expression/genes && \
echo -e "align_and_estimate_abundance/GM12878_Rep1/GM12878_Rep1.genes.results\nalign_and_estimate_abundance/GM12878_Rep2/GM12878_Rep2.genes.results\nalign_and_estimate_abundance/H1ESC_Rep1/H1ESC_Rep1.genes.results\nalign_and_estimate_abundance/H1ESC_Rep2/H1ESC_Rep2.genes.results" > differential_expression/genes.counts.files && \
abundance_estimates_to_matrix.pl \
  --est_method RSEM \
  --out_prefix differential_expression/genes \
  differential_expression/genes.counts.files && \
awk -F '\t' '{OFS="\t" ; print $1,$0}' differential_expression/genes.counts.matrix | sed '1s/^\t/Genes\tSymbol/' \
  > differential_expression/genes.counts.matrix.symbol && \
cut -f 1,3,4 align_and_estimate_abundance/GM12878_Rep1/GM12878_Rep1.genes.results > differential_expression/genes.lengths.tsv &&  sed '1d' differential_expression/genes.lengths.tsv > differential_expression/genes.lengths.tsv.noheader.tsv
align_and_estimate_abundance.genes.b3c4c1fd6af7e4eaf50b3947423ed48a.mugqic.done
chmod 755 $COMMAND
align_and_estimate_abundance_5_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&     $COMMAND 
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE


if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem-per-cpu 4700M -c 5 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$align_and_estimate_abundance_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$align_and_estimate_abundance_5_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: align_and_estimate_abundance_6_JOB_ID: align_and_estimate_abundance.isoforms
#-------------------------------------------------------------------------------
JOB_NAME=align_and_estimate_abundance.isoforms
JOB_DEPENDENCIES=$align_and_estimate_abundance_1_JOB_ID:$align_and_estimate_abundance_2_JOB_ID:$align_and_estimate_abundance_3_JOB_ID:$align_and_estimate_abundance_4_JOB_ID
JOB_DONE=job_output/align_and_estimate_abundance/align_and_estimate_abundance.isoforms.5490fd0ab6e6e2aeee750c0cfa5b1f9e.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'align_and_estimate_abundance.isoforms.5490fd0ab6e6e2aeee750c0cfa5b1f9e.mugqic.done' > $COMMAND
module purge && \
module load mugqic/perl/5.22.1 mugqic/trinity/2.0.4_patch mugqic/R_Bioconductor/3.6.0_3.9 && \
mkdir -p differential_expression/isoforms && \
echo -e "align_and_estimate_abundance/GM12878_Rep1/GM12878_Rep1.isoforms.results\nalign_and_estimate_abundance/GM12878_Rep2/GM12878_Rep2.isoforms.results\nalign_and_estimate_abundance/H1ESC_Rep1/H1ESC_Rep1.isoforms.results\nalign_and_estimate_abundance/H1ESC_Rep2/H1ESC_Rep2.isoforms.results" > differential_expression/isoforms.counts.files && \
abundance_estimates_to_matrix.pl \
  --est_method RSEM \
  --out_prefix differential_expression/isoforms \
  differential_expression/isoforms.counts.files && \
awk -F '\t' '{OFS="\t" ; print $1,$0}' differential_expression/isoforms.counts.matrix | sed '1s/^\t/Isoforms\tSymbol/' \
  > differential_expression/isoforms.counts.matrix.symbol && \
cut -f 1,3,4 align_and_estimate_abundance/GM12878_Rep1/GM12878_Rep1.isoforms.results > differential_expression/isoforms.lengths.tsv &&  sed '1d' differential_expression/isoforms.lengths.tsv > differential_expression/isoforms.lengths.tsv.noheader.tsv
align_and_estimate_abundance.isoforms.5490fd0ab6e6e2aeee750c0cfa5b1f9e.mugqic.done
chmod 755 $COMMAND
align_and_estimate_abundance_6_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&     $COMMAND 
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE


if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem-per-cpu 4700M -c 5 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$align_and_estimate_abundance_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$align_and_estimate_abundance_6_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: align_and_estimate_abundance_7_JOB_ID: align_and_estimate_abundance.parse_trinotate
#-------------------------------------------------------------------------------
JOB_NAME=align_and_estimate_abundance.parse_trinotate
JOB_DEPENDENCIES=$trinotate_1_JOB_ID:$align_and_estimate_abundance_6_JOB_ID
JOB_DONE=job_output/align_and_estimate_abundance/align_and_estimate_abundance.parse_trinotate.3402a69b523b1aaaf48d365a4b784a43.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'align_and_estimate_abundance.parse_trinotate.3402a69b523b1aaaf48d365a4b784a43.mugqic.done' > $COMMAND
module purge && \
module load mugqic/mugqic_tools/2.8.2 mugqic/python/2.7.13 && \
$PYTHON_TOOLS/parseTrinotateOutput.py -r trinotate/trinotate_annotation_report.tsv -o trinotate/trinotate_annotation_report.tsv.genes -i "#gene_id" -l differential_expression/isoforms.lengths.tsv &&
$PYTHON_TOOLS/parseTrinotateOutput.py -r trinotate/trinotate_annotation_report.tsv -o trinotate/trinotate_annotation_report.tsv.isoforms -i "transcript_id" -f sprot_Top_BLASTX_hit != \".\" or TmHMM != \".\"
align_and_estimate_abundance.parse_trinotate.3402a69b523b1aaaf48d365a4b784a43.mugqic.done
chmod 755 $COMMAND
align_and_estimate_abundance_7_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&     $COMMAND 
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE


if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem-per-cpu 4700M -c 5 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$align_and_estimate_abundance_7_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$align_and_estimate_abundance_7_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# STEP: gq_seq_utils_exploratory_analysis_rnaseq_denovo
#-------------------------------------------------------------------------------
STEP=gq_seq_utils_exploratory_analysis_rnaseq_denovo
mkdir -p $JOB_OUTPUT_DIR/$STEP

#-------------------------------------------------------------------------------
# JOB: gq_seq_utils_exploratory_analysis_rnaseq_denovo_1_JOB_ID: gq_seq_utils_exploratory_analysis_rnaseq_denovo
#-------------------------------------------------------------------------------
JOB_NAME=gq_seq_utils_exploratory_analysis_rnaseq_denovo
JOB_DEPENDENCIES=$align_and_estimate_abundance_5_JOB_ID
JOB_DONE=job_output/gq_seq_utils_exploratory_analysis_rnaseq_denovo/gq_seq_utils_exploratory_analysis_rnaseq_denovo.bfbd4de4b11b8e0063c9563c864a3092.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gq_seq_utils_exploratory_analysis_rnaseq_denovo.bfbd4de4b11b8e0063c9563c864a3092.mugqic.done' > $COMMAND
module purge && \
module load mugqic/R_Bioconductor/3.6.0_3.9 mugqic/mugqic_R_packages/1.0.5 && \
mkdir -p exploratory && \
R --no-save --no-restore <<-EOF
suppressPackageStartupMessages(library(gqSeqUtils))
exploratoryAnalysisRNAseqdenovo(read.counts.path="differential_expression/genes.counts.matrix", genes.path="differential_expression/genes.lengths.tsv", output.dir="exploratory")
desc = readRDS(file.path("exploratory","index.RData"))
write.table(desc,file=file.path("exploratory","index.tsv"),sep='\t',quote=F,col.names=T,row.names=F)
print("done.")

EOF
gq_seq_utils_exploratory_analysis_rnaseq_denovo.bfbd4de4b11b8e0063c9563c864a3092.mugqic.done
chmod 755 $COMMAND
gq_seq_utils_exploratory_analysis_rnaseq_denovo_1_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&     $COMMAND 
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE


if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=00:30:00 --mem-per-cpu 4700M -c 1 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$gq_seq_utils_exploratory_analysis_rnaseq_denovo_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$gq_seq_utils_exploratory_analysis_rnaseq_denovo_1_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: gq_seq_utils_exploratory_analysis_rnaseq_denovo_2_JOB_ID: gq_seq_utils_exploratory_analysis_rnaseq_denovo_report
#-------------------------------------------------------------------------------
JOB_NAME=gq_seq_utils_exploratory_analysis_rnaseq_denovo_report
JOB_DEPENDENCIES=$gq_seq_utils_exploratory_analysis_rnaseq_denovo_1_JOB_ID
JOB_DONE=job_output/gq_seq_utils_exploratory_analysis_rnaseq_denovo/gq_seq_utils_exploratory_analysis_rnaseq_denovo_report.23219cfe2040024f5cf0d88d305f75c0.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gq_seq_utils_exploratory_analysis_rnaseq_denovo_report.23219cfe2040024f5cf0d88d305f75c0.mugqic.done' > $COMMAND
module purge && \
module load mugqic/R_Bioconductor/3.6.0_3.9 mugqic/pandoc/2.16.1 && \
R --no-save --no-restore <<-'EOF'
report_dir="report";
input_rmarkdown_file = '/cvmfs/soft.mugqic/root/software/genpipes/genpipes-4.3.1/bfx/report/RnaSeqDeNovoAssembly.gq_seq_utils_exploratory_analysis_rnaseq.Rmd'
render_output_dir    = 'report'
rmarkdown_file       = basename(input_rmarkdown_file) # honoring a different WD that location of Rmd file in knitr is problematic
file.copy(from = input_rmarkdown_file, to = rmarkdown_file, overwrite = T)
rmarkdown::render(input = rmarkdown_file, output_format = c("html_document","md_document"), output_dir = render_output_dir  )
file.remove(rmarkdown_file)
EOF
gq_seq_utils_exploratory_analysis_rnaseq_denovo_report.23219cfe2040024f5cf0d88d305f75c0.mugqic.done
chmod 755 $COMMAND
gq_seq_utils_exploratory_analysis_rnaseq_denovo_2_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&     $COMMAND 
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE


if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem-per-cpu 4700M -c 1 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$gq_seq_utils_exploratory_analysis_rnaseq_denovo_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$gq_seq_utils_exploratory_analysis_rnaseq_denovo_2_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# STEP: differential_expression
#-------------------------------------------------------------------------------
STEP=differential_expression
mkdir -p $JOB_OUTPUT_DIR/$STEP

#-------------------------------------------------------------------------------
# JOB: differential_expression_1_JOB_ID: differential_expression_genes
#-------------------------------------------------------------------------------
JOB_NAME=differential_expression_genes
JOB_DEPENDENCIES=$trinotate_1_JOB_ID:$align_and_estimate_abundance_5_JOB_ID:$align_and_estimate_abundance_7_JOB_ID
JOB_DONE=job_output/differential_expression/differential_expression_genes.76bb4de9b06f67ff725905d7c7f77a4d.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'differential_expression_genes.76bb4de9b06f67ff725905d7c7f77a4d.mugqic.done' > $COMMAND
module purge && \
module load mugqic/mugqic_tools/2.8.2 mugqic/R_Bioconductor/3.6.0_3.9 mugqic/python/2.7.13 && \
mkdir -p differential_expression/genes && \
Rscript $R_TOOLS/edger.R \
  -d design.rnaseq.txt \
  -c differential_expression/genes.counts.matrix.symbol \
  -o differential_expression/genes && \
Rscript $R_TOOLS/deseq2.R \
  -d design.rnaseq.txt \
  -c differential_expression/genes.counts.matrix.symbol \
  -o differential_expression/genes \
   && \
$PYTHON_TOOLS/parseMergeCsv.py -i differential_expression/genes/H1ESC_GM12787/dge_results.csv trinotate/trinotate_annotation_report.tsv.genes_blast.tsv \
      -o differential_expression/genes/H1ESC_GM12787/dge_trinotate_results.csv \
      -c id "#gene_id" \
      -d \\t  -x "#gene_id" transcript_id -l  -t edger.p.value -n  && \
Rscript $R_TOOLS/goseq.R -p 0.1 -f 0.1 \
  -a differential_expression/genes.lengths.tsv.noheader.tsv \
  -G trinotate/trinotate_annotation_report.tsv.genes_go.tsv \
  -d differential_expression/genes/H1ESC_GM12787/dge_trinotate_results.csv \
  -c 1,6 \
  -o differential_expression/genes/H1ESC_GM12787/gene_ontology_results.csv
differential_expression_genes.76bb4de9b06f67ff725905d7c7f77a4d.mugqic.done
chmod 755 $COMMAND
differential_expression_1_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&     $COMMAND 
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE


if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem-per-cpu 4700M -c 1 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$differential_expression_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$differential_expression_1_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: differential_expression_2_JOB_ID: differential_expression_isoforms
#-------------------------------------------------------------------------------
JOB_NAME=differential_expression_isoforms
JOB_DEPENDENCIES=$trinotate_1_JOB_ID:$align_and_estimate_abundance_6_JOB_ID:$align_and_estimate_abundance_7_JOB_ID
JOB_DONE=job_output/differential_expression/differential_expression_isoforms.c185887bfcf93570570bf084e9625f62.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'differential_expression_isoforms.c185887bfcf93570570bf084e9625f62.mugqic.done' > $COMMAND
module purge && \
module load mugqic/mugqic_tools/2.8.2 mugqic/R_Bioconductor/3.6.0_3.9 mugqic/python/2.7.13 && \
mkdir -p differential_expression/isoforms && \
Rscript $R_TOOLS/edger.R \
  -d design.rnaseq.txt \
  -c differential_expression/isoforms.counts.matrix.symbol \
  -o differential_expression/isoforms && \
Rscript $R_TOOLS/deseq2.R \
  -d design.rnaseq.txt \
  -c differential_expression/isoforms.counts.matrix.symbol \
  -o differential_expression/isoforms \
   && \
$PYTHON_TOOLS/parseMergeCsv.py -i differential_expression/isoforms/H1ESC_GM12787/dge_results.csv trinotate/trinotate_annotation_report.tsv.isoforms_blast.tsv \
      -o differential_expression/isoforms/H1ESC_GM12787/dge_trinotate_results.csv \
      -c id transcript_id \
      -d \\t  -x "#gene_id" transcript_id -l  -t edger.p.value -n  && \
Rscript $R_TOOLS/goseq.R -p 0.1 -f 0.1 \
  -a differential_expression/isoforms.lengths.tsv.noheader.tsv \
  -G trinotate/trinotate_annotation_report.tsv.isoforms_go.tsv \
  -d differential_expression/isoforms/H1ESC_GM12787/dge_trinotate_results.csv \
  -c 1,6 \
  -o differential_expression/isoforms/H1ESC_GM12787/gene_ontology_results.csv
differential_expression_isoforms.c185887bfcf93570570bf084e9625f62.mugqic.done
chmod 755 $COMMAND
differential_expression_2_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&     $COMMAND 
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE


if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem-per-cpu 4700M -c 1 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$differential_expression_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$differential_expression_2_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: differential_expression_3_JOB_ID: differential_expression_goseq_rnaseq_denovo_report
#-------------------------------------------------------------------------------
JOB_NAME=differential_expression_goseq_rnaseq_denovo_report
JOB_DEPENDENCIES=$differential_expression_1_JOB_ID:$differential_expression_2_JOB_ID
JOB_DONE=job_output/differential_expression/differential_expression_goseq_rnaseq_denovo_report.9dd93e8ec295b9d8dbb492493126bd83.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'differential_expression_goseq_rnaseq_denovo_report.9dd93e8ec295b9d8dbb492493126bd83.mugqic.done' > $COMMAND
module purge && \
module load mugqic/R_Bioconductor/3.6.0_3.9 mugqic/pandoc/2.16.1 && \
R --no-save --no-restore <<-'EOF'
design_file="design.rnaseq.txt"; report_dir="report"; source_dir="differential_expression"; top_n_results=10; contrasts=c("H1ESC_GM12787");
input_rmarkdown_file = '/cvmfs/soft.mugqic/root/software/genpipes/genpipes-4.3.1/bfx/report/RnaSeqDeNovoAssembly.differential_expression_goseq.Rmd'
render_output_dir    = 'report'
rmarkdown_file       = basename(input_rmarkdown_file) # honoring a different WD that location of Rmd file in knitr is problematic
file.copy(from = input_rmarkdown_file, to = rmarkdown_file, overwrite = T)
rmarkdown::render(input = rmarkdown_file, output_format = c("html_document","md_document"), output_dir = render_output_dir  )
file.remove(rmarkdown_file)
EOF
differential_expression_goseq_rnaseq_denovo_report.9dd93e8ec295b9d8dbb492493126bd83.mugqic.done
chmod 755 $COMMAND
differential_expression_3_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&     $COMMAND 
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE


if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem-per-cpu 4700M -c 1 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$differential_expression_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$differential_expression_3_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# STEP: filter_annotated_components
#-------------------------------------------------------------------------------
STEP=filter_annotated_components
mkdir -p $JOB_OUTPUT_DIR/$STEP

#-------------------------------------------------------------------------------
# JOB: filter_annotated_components_1_JOB_ID: filter_annotated_components
#-------------------------------------------------------------------------------
JOB_NAME=filter_annotated_components
JOB_DEPENDENCIES=$trinity_1_JOB_ID:$align_and_estimate_abundance_7_JOB_ID
JOB_DONE=job_output/filter_annotated_components/filter_annotated_components.1fb7e7b014ed9fcf7c8ac27d12b563a2.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'filter_annotated_components.1fb7e7b014ed9fcf7c8ac27d12b563a2.mugqic.done' > $COMMAND
module purge && \
module load mugqic/mugqic_tools/2.8.2 mugqic/python/2.7.13 mugqic/R_Bioconductor/3.6.0_3.9 mugqic/mugqic_R_packages/1.0.6 && \
mkdir -p filtered_assembly && \
$PYTHON_TOOLS/filterAssemblyToFastaToXls.py -f trinity_out_dir/Trinity.fasta \
-o filtered_assembly/Trinity \
-l trinotate/trinotate_annotation_report.tsv.isoforms_filtered.tsv \
-c 0  && \
Rscript -e 'library(gqSeqUtils);dnaFastaStats(filename="filtered_assembly/Trinity.fasta",type="trinity",output.prefix="filtered_assembly/trinity_filtered.stats")' && \
zip -j filtered_assembly/Trinity.fasta.zip filtered_assembly/Trinity.fasta filtered_assembly/Trinity.tsv
filter_annotated_components.1fb7e7b014ed9fcf7c8ac27d12b563a2.mugqic.done
chmod 755 $COMMAND
filter_annotated_components_1_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&     $COMMAND 
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE


if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem-per-cpu 4700M -c 1 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$filter_annotated_components_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$filter_annotated_components_1_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: filter_annotated_components_2_JOB_ID: filter_annotated_components_report
#-------------------------------------------------------------------------------
JOB_NAME=filter_annotated_components_report
JOB_DEPENDENCIES=$filter_annotated_components_1_JOB_ID
JOB_DONE=job_output/filter_annotated_components/filter_annotated_components_report.3967e253de030ee54c0a0c19ca79fbce.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'filter_annotated_components_report.3967e253de030ee54c0a0c19ca79fbce.mugqic.done' > $COMMAND
module purge && \
module load mugqic/pandoc/2.16.1 && \
mkdir -p report && \
cp filtered_assembly/Trinity.fasta.zip report/filtered_assembly.zip && \
cp filtered_assembly/trinity_filtered.stats.csv filtered_assembly/trinity_filtered.stats.jpg filtered_assembly/trinity_filtered.stats.pdf report/ && \
assembly_table=`sed '1d' filtered_assembly/trinity_filtered.stats.csv | perl -pe 's/^"([^"]*)",/\1	/g' | grep -P "^(Nb. Transcripts|Nb. Components|Total Transcripts Length|Min. Transcript Length|Median Transcript Length|Mean Transcript Length|Max. Transcript Length|N50)" | LC_NUMERIC=en_CA awk -F"	" '{print $1"|"sprintf("%\47d", $2)}'` && \
pandoc --to=markdown \
--template /cvmfs/soft.mugqic/root/software/genpipes/genpipes-4.3.1/bfx/report/RnaSeqDeNovoAssembly.filtered.trinity.md \
--variable assembly_table="$assembly_table" \
--variable filter_string="sprot_Top_BLASTX_hit != \".\" or TmHMM != \".\"" \
/cvmfs/soft.mugqic/root/software/genpipes/genpipes-4.3.1/bfx/report/RnaSeqDeNovoAssembly.filtered.trinity.md \
> report/RnaSeqDeNovoAssembly.filtered.trinity.md
filter_annotated_components_report.3967e253de030ee54c0a0c19ca79fbce.mugqic.done
chmod 755 $COMMAND
filter_annotated_components_2_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&     $COMMAND 
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE


if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem-per-cpu 4700M -c 1 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$filter_annotated_components_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$filter_annotated_components_2_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# STEP: gq_seq_utils_exploratory_analysis_rnaseq_denovo_filtered
#-------------------------------------------------------------------------------
STEP=gq_seq_utils_exploratory_analysis_rnaseq_denovo_filtered
mkdir -p $JOB_OUTPUT_DIR/$STEP

#-------------------------------------------------------------------------------
# JOB: gq_seq_utils_exploratory_analysis_rnaseq_denovo_filtered_1_JOB_ID: filter_annotated_components_exploratory
#-------------------------------------------------------------------------------
JOB_NAME=filter_annotated_components_exploratory
JOB_DEPENDENCIES=$align_and_estimate_abundance_6_JOB_ID:$align_and_estimate_abundance_7_JOB_ID
JOB_DONE=job_output/gq_seq_utils_exploratory_analysis_rnaseq_denovo_filtered/filter_annotated_components_exploratory.3a666df90f22b59436d5c477dea24643.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'filter_annotated_components_exploratory.3a666df90f22b59436d5c477dea24643.mugqic.done' > $COMMAND
module purge && \
module load mugqic/mugqic_tools/2.8.2 mugqic/python/2.7.13 && \
mkdir -p filtered_assembly/exploratory && \
sed '1s/^/ \n/' trinotate/trinotate_annotation_report.tsv.isoforms_filtered.tsv > trinotate/trinotate_annotation_report.tsv.isoforms_filtered_header.tsv && \
$PYTHON_TOOLS/parseMergeCsv.py -i trinotate/trinotate_annotation_report.tsv.isoforms_filtered_header.tsv differential_expression/isoforms.counts.matrix \
      -o filtered_assembly/isoforms.counts.matrix \
      -c '' \
      -d \\t  -x '' -l  && \
$PYTHON_TOOLS/parseMergeCsv.py -i trinotate/trinotate_annotation_report.tsv.isoforms_filtered_header.tsv differential_expression/isoforms.lengths.tsv \
      -o filtered_assembly/isoforms.lengths.tsv \
      -c '' transcript_id \
      -d \\t  -x ' ' -l 
filter_annotated_components_exploratory.3a666df90f22b59436d5c477dea24643.mugqic.done
chmod 755 $COMMAND
gq_seq_utils_exploratory_analysis_rnaseq_denovo_filtered_1_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&     $COMMAND 
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE


if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem-per-cpu 4700M -c 1 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$gq_seq_utils_exploratory_analysis_rnaseq_denovo_filtered_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$gq_seq_utils_exploratory_analysis_rnaseq_denovo_filtered_1_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: gq_seq_utils_exploratory_analysis_rnaseq_denovo_filtered_2_JOB_ID: gq_seq_utils_exploratory_analysis_rnaseq_denovo
#-------------------------------------------------------------------------------
JOB_NAME=gq_seq_utils_exploratory_analysis_rnaseq_denovo
JOB_DEPENDENCIES=$gq_seq_utils_exploratory_analysis_rnaseq_denovo_filtered_1_JOB_ID
JOB_DONE=job_output/gq_seq_utils_exploratory_analysis_rnaseq_denovo_filtered/gq_seq_utils_exploratory_analysis_rnaseq_denovo.cff6cc971c99b7eee6a04ca0010a5d73.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gq_seq_utils_exploratory_analysis_rnaseq_denovo.cff6cc971c99b7eee6a04ca0010a5d73.mugqic.done' > $COMMAND
module purge && \
module load mugqic/R_Bioconductor/3.6.0_3.9 mugqic/mugqic_R_packages/1.0.5 && \
mkdir -p filtered_assembly/exploratory && \
R --no-save --no-restore <<-EOF
suppressPackageStartupMessages(library(gqSeqUtils))
exploratoryAnalysisRNAseqdenovo(read.counts.path="filtered_assembly/isoforms.counts.matrix", genes.path="filtered_assembly/isoforms.lengths.tsv", output.dir="filtered_assembly/exploratory")
desc = readRDS(file.path("filtered_assembly/exploratory","index.RData"))
write.table(desc,file=file.path("filtered_assembly/exploratory","index.tsv"),sep='\t',quote=F,col.names=T,row.names=F)
print("done.")

EOF
gq_seq_utils_exploratory_analysis_rnaseq_denovo.cff6cc971c99b7eee6a04ca0010a5d73.mugqic.done
chmod 755 $COMMAND
gq_seq_utils_exploratory_analysis_rnaseq_denovo_filtered_2_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&     $COMMAND 
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE


if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=00:30:00 --mem-per-cpu 4700M -c 1 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$gq_seq_utils_exploratory_analysis_rnaseq_denovo_filtered_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$gq_seq_utils_exploratory_analysis_rnaseq_denovo_filtered_2_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: gq_seq_utils_exploratory_analysis_rnaseq_denovo_filtered_3_JOB_ID: gq_seq_utils_exploratory_analysis_rnaseq_denovo_filtered_report
#-------------------------------------------------------------------------------
JOB_NAME=gq_seq_utils_exploratory_analysis_rnaseq_denovo_filtered_report
JOB_DEPENDENCIES=$gq_seq_utils_exploratory_analysis_rnaseq_denovo_filtered_2_JOB_ID
JOB_DONE=job_output/gq_seq_utils_exploratory_analysis_rnaseq_denovo_filtered/gq_seq_utils_exploratory_analysis_rnaseq_denovo_filtered_report.106682e3d99013eaf7a3187c205a32f6.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gq_seq_utils_exploratory_analysis_rnaseq_denovo_filtered_report.106682e3d99013eaf7a3187c205a32f6.mugqic.done' > $COMMAND
module purge && \
module load mugqic/R_Bioconductor/3.6.0_3.9 mugqic/pandoc/2.16.1 && \
R --no-save --no-restore <<-'EOF'
report_dir="report/filtered_assembly"; exploratory_dir="filtered_assembly/exploratory";
input_rmarkdown_file = '/cvmfs/soft.mugqic/root/software/genpipes/genpipes-4.3.1/bfx/report/RnaSeqDeNovoAssembly.gq_seq_utils_exploratory_analysis_rnaseq_filtered.Rmd'
render_output_dir    = 'report'
rmarkdown_file       = basename(input_rmarkdown_file) # honoring a different WD that location of Rmd file in knitr is problematic
file.copy(from = input_rmarkdown_file, to = rmarkdown_file, overwrite = T)
rmarkdown::render(input = rmarkdown_file, output_format = c("html_document","md_document"), output_dir = render_output_dir  )
file.remove(rmarkdown_file)
EOF
gq_seq_utils_exploratory_analysis_rnaseq_denovo_filtered_report.106682e3d99013eaf7a3187c205a32f6.mugqic.done
chmod 755 $COMMAND
gq_seq_utils_exploratory_analysis_rnaseq_denovo_filtered_3_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&     $COMMAND 
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE


if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem-per-cpu 4700M -c 1 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$gq_seq_utils_exploratory_analysis_rnaseq_denovo_filtered_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$gq_seq_utils_exploratory_analysis_rnaseq_denovo_filtered_3_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# STEP: differential_expression_filtered
#-------------------------------------------------------------------------------
STEP=differential_expression_filtered
mkdir -p $JOB_OUTPUT_DIR/$STEP

#-------------------------------------------------------------------------------
# JOB: differential_expression_filtered_1_JOB_ID: differential_expression_filtered_get_trinotate
#-------------------------------------------------------------------------------
JOB_NAME=differential_expression_filtered_get_trinotate
JOB_DEPENDENCIES=$align_and_estimate_abundance_7_JOB_ID
JOB_DONE=job_output/differential_expression_filtered/differential_expression_filtered_get_trinotate.363b87286693b6df3d6c9f87945f6bda.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'differential_expression_filtered_get_trinotate.363b87286693b6df3d6c9f87945f6bda.mugqic.done' > $COMMAND
mkdir -p filtered_assembly/differential_expression && \
cat trinotate/trinotate_annotation_report.tsv.isoforms_filtered.tsv | awk 'BEGIN{OFS="_";FS="_"}{print $1,$2}' | uniq | sed '1s/^/ \n/'   > trinotate/trinotate_annotation_report.tsv.genes_filtered_header.tsv && \
sed '1s/^/ \n/' trinotate/trinotate_annotation_report.tsv.isoforms_filtered.tsv > trinotate/trinotate_annotation_report.tsv.isoforms_filtered_header.tsv
differential_expression_filtered_get_trinotate.363b87286693b6df3d6c9f87945f6bda.mugqic.done
chmod 755 $COMMAND
differential_expression_filtered_1_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&     $COMMAND 
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE


if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem-per-cpu 4700M -c 1 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$differential_expression_filtered_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$differential_expression_filtered_1_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: differential_expression_filtered_2_JOB_ID: differential_expression_filtered_genes
#-------------------------------------------------------------------------------
JOB_NAME=differential_expression_filtered_genes
JOB_DEPENDENCIES=$trinotate_1_JOB_ID:$align_and_estimate_abundance_5_JOB_ID:$align_and_estimate_abundance_7_JOB_ID:$differential_expression_filtered_1_JOB_ID
JOB_DONE=job_output/differential_expression_filtered/differential_expression_filtered_genes.c7fbfd0ce19d1ca7a7abe23a732e7480.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'differential_expression_filtered_genes.c7fbfd0ce19d1ca7a7abe23a732e7480.mugqic.done' > $COMMAND
module purge && \
module load mugqic/mugqic_tools/2.8.2 mugqic/python/2.7.13 mugqic/R_Bioconductor/3.6.0_3.9 && \
$PYTHON_TOOLS/parseMergeCsv.py -i trinotate/trinotate_annotation_report.tsv.genes_filtered_header.tsv differential_expression/genes.counts.matrix.symbol \
      -o filtered_assembly/differential_expression/genes.counts.matrix.symbol \
      -c '' Genes \
      -d \\t  -x ' ' -l  && \
cp differential_expression/genes.lengths.tsv.noheader.tsv filtered_assembly/differential_expression/genes.lengths.tsv.noheader.tsv && \
mkdir -p filtered_assembly/differential_expression/genes && \
Rscript $R_TOOLS/edger.R \
  -d design.rnaseq.txt \
  -c filtered_assembly/differential_expression/genes.counts.matrix.symbol \
  -o filtered_assembly/differential_expression/genes && \
Rscript $R_TOOLS/deseq2.R \
  -d design.rnaseq.txt \
  -c filtered_assembly/differential_expression/genes.counts.matrix.symbol \
  -o filtered_assembly/differential_expression/genes \
   && \
$PYTHON_TOOLS/parseMergeCsv.py -i filtered_assembly/differential_expression/genes/H1ESC_GM12787/dge_results.csv trinotate/trinotate_annotation_report.tsv.genes_blast.tsv \
      -o filtered_assembly/differential_expression/genes/H1ESC_GM12787/dge_trinotate_results.csv \
      -c id "#gene_id" \
      -d \\t  -x "#gene_id" transcript_id -l  -t edger.p.value -n  && \
Rscript $R_TOOLS/goseq.R -p 0.1 -f 0.1 \
  -a filtered_assembly/differential_expression/genes.lengths.tsv.noheader.tsv \
  -G trinotate/trinotate_annotation_report.tsv.genes_go.tsv \
  -d filtered_assembly/differential_expression/genes/H1ESC_GM12787/dge_trinotate_results.csv \
  -c 1,6 \
  -o filtered_assembly/differential_expression/genes/H1ESC_GM12787/gene_ontology_results.csv
differential_expression_filtered_genes.c7fbfd0ce19d1ca7a7abe23a732e7480.mugqic.done
chmod 755 $COMMAND
differential_expression_filtered_2_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&     $COMMAND 
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE


if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem-per-cpu 4700M -c 1 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$differential_expression_filtered_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$differential_expression_filtered_2_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: differential_expression_filtered_3_JOB_ID: differential_expression_filtered_isoforms
#-------------------------------------------------------------------------------
JOB_NAME=differential_expression_filtered_isoforms
JOB_DEPENDENCIES=$trinotate_1_JOB_ID:$align_and_estimate_abundance_6_JOB_ID:$align_and_estimate_abundance_7_JOB_ID:$gq_seq_utils_exploratory_analysis_rnaseq_denovo_filtered_1_JOB_ID:$differential_expression_filtered_1_JOB_ID
JOB_DONE=job_output/differential_expression_filtered/differential_expression_filtered_isoforms.d9512045dec0dddbcd2487710b2a3e19.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'differential_expression_filtered_isoforms.d9512045dec0dddbcd2487710b2a3e19.mugqic.done' > $COMMAND
module purge && \
module load mugqic/mugqic_tools/2.8.2 mugqic/python/2.7.13 mugqic/R_Bioconductor/3.6.0_3.9 && \
$PYTHON_TOOLS/parseMergeCsv.py -i trinotate/trinotate_annotation_report.tsv.isoforms_filtered_header.tsv differential_expression/isoforms.counts.matrix.symbol \
      -o filtered_assembly/differential_expression/isoforms.counts.matrix.symbol \
      -c '' Isoforms \
      -d \\t  -x ' ' -l  && \
cp differential_expression/isoforms.lengths.tsv.noheader.tsv filtered_assembly/differential_expression/isoforms.lengths.tsv.noheader.tsv && \
mkdir -p filtered_assembly/differential_expression/isoforms && \
Rscript $R_TOOLS/edger.R \
  -d design.rnaseq.txt \
  -c filtered_assembly/differential_expression/isoforms.counts.matrix.symbol \
  -o filtered_assembly/differential_expression/isoforms && \
Rscript $R_TOOLS/deseq2.R \
  -d design.rnaseq.txt \
  -c filtered_assembly/differential_expression/isoforms.counts.matrix.symbol \
  -o filtered_assembly/differential_expression/isoforms \
   && \
$PYTHON_TOOLS/parseMergeCsv.py -i filtered_assembly/differential_expression/isoforms/H1ESC_GM12787/dge_results.csv trinotate/trinotate_annotation_report.tsv.isoforms_blast.tsv \
      -o filtered_assembly/differential_expression/isoforms/H1ESC_GM12787/dge_trinotate_results.csv \
      -c id transcript_id \
      -d \\t  -x "#gene_id" transcript_id -l  -t edger.p.value -n  && \
Rscript $R_TOOLS/goseq.R -p 0.1 -f 0.1 \
  -a filtered_assembly/differential_expression/isoforms.lengths.tsv.noheader.tsv \
  -G trinotate/trinotate_annotation_report.tsv.isoforms_go.tsv \
  -d filtered_assembly/differential_expression/isoforms/H1ESC_GM12787/dge_trinotate_results.csv \
  -c 1,6 \
  -o filtered_assembly/differential_expression/isoforms/H1ESC_GM12787/gene_ontology_results.csv
differential_expression_filtered_isoforms.d9512045dec0dddbcd2487710b2a3e19.mugqic.done
chmod 755 $COMMAND
differential_expression_filtered_3_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&     $COMMAND 
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE


if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem-per-cpu 4700M -c 1 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$differential_expression_filtered_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$differential_expression_filtered_3_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: differential_expression_filtered_4_JOB_ID: differential_expression_goseq_rnaseq_denovo_filtered_report
#-------------------------------------------------------------------------------
JOB_NAME=differential_expression_goseq_rnaseq_denovo_filtered_report
JOB_DEPENDENCIES=$gq_seq_utils_exploratory_analysis_rnaseq_denovo_filtered_1_JOB_ID:$differential_expression_filtered_1_JOB_ID:$differential_expression_filtered_2_JOB_ID:$differential_expression_filtered_3_JOB_ID
JOB_DONE=job_output/differential_expression_filtered/differential_expression_goseq_rnaseq_denovo_filtered_report.921b35317ca0a6bb32128bd5bd6d27c8.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'differential_expression_goseq_rnaseq_denovo_filtered_report.921b35317ca0a6bb32128bd5bd6d27c8.mugqic.done' > $COMMAND
module purge && \
module load mugqic/R_Bioconductor/3.6.0_3.9 mugqic/pandoc/2.16.1 && \
R --no-save --no-restore <<-'EOF'
report_dir="report/filtered_assembly"; source_dir="filtered_assembly/differential_expression"; top_n_results=10; contrasts=c("H1ESC_GM12787");
input_rmarkdown_file = '/cvmfs/soft.mugqic/root/software/genpipes/genpipes-4.3.1/bfx/report/RnaSeqDeNovoAssembly.differential_expression_goseq_filtered.Rmd'
render_output_dir    = 'report'
rmarkdown_file       = basename(input_rmarkdown_file) # honoring a different WD that location of Rmd file in knitr is problematic
file.copy(from = input_rmarkdown_file, to = rmarkdown_file, overwrite = T)
rmarkdown::render(input = rmarkdown_file, output_format = c("html_document","md_document"), output_dir = render_output_dir  )
file.remove(rmarkdown_file)
EOF
differential_expression_goseq_rnaseq_denovo_filtered_report.921b35317ca0a6bb32128bd5bd6d27c8.mugqic.done
chmod 755 $COMMAND
differential_expression_filtered_4_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&     $COMMAND 
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE


if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem-per-cpu 4700M -c 1 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$differential_expression_filtered_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$differential_expression_filtered_4_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# Call home with pipeline statistics
#-------------------------------------------------------------------------------
LOG_MD5=$(echo $USER-'10.74.73.1-RnaSeqDeNovoAssembly-GM12878_Rep1.GM12878_Rep1,GM12878_Rep2.GM12878_Rep2,H1ESC_Rep1.H1ESC_Rep1,H1ESC_Rep2.H1ESC_Rep2' | md5sum | awk '{ print $1 }')
if test -t 1; then ncolors=$(tput colors); if test -n "$ncolors" && test $ncolors -ge 8; then bold="$(tput bold)"; normal="$(tput sgr0)"; yellow="$(tput setaf 3)"; fi; fi
wget --quiet 'http://mugqic.hpc.mcgill.ca/cgi-bin/pipeline.cgi?hostname=beluga1.int.ets1.calculquebec.ca&ip=10.74.73.1&pipeline=RnaSeqDeNovoAssembly&steps=picard_sam_to_fastq,trimmomatic,merge_trimmomatic_stats,insilico_read_normalization_readsets,insilico_read_normalization_all,trinity,exonerate_fastasplit,blastx_trinity_uniprot,blastx_trinity_uniprot_merge,transdecoder,hmmer,rnammer_transcriptome,blastp_transdecoder_uniprot,signalp,tmhmm,trinotate,align_and_estimate_abundance_prep_reference,align_and_estimate_abundance,gq_seq_utils_exploratory_analysis_rnaseq_denovo,differential_expression,filter_annotated_components,gq_seq_utils_exploratory_analysis_rnaseq_denovo_filtered,differential_expression_filtered&samples=4&md5=$LOG_MD5' -O /dev/null || echo "${bold}${yellow}Warning:${normal}${yellow} Genpipes ran successfully but was not send telemetry to mugqic.hpc.mcgill.ca. This error will not affect genpipes jobs you have submitted.${normal}"
