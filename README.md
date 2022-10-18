## __denovo_RNAseq_GenPipes_Workshop2022__

___

 Demo data , slides and protocol to run Trinity based *de novo* RNA (transcriptome) assembly pipeline using GenPipes. 

___

Author: Senthilkumar Kailasam

Date: 17 Oct, 2022

<br>
<br>

### __1. Install all  prerequisites/dependency__

We  will go through the steps involved in installing Genpipes in Ubuntu 22.04 LTS


Tools required:

- tar 
- podman (https://podman.io/getting-started/installation)
- git (https://docs.gitlab.com/ee/topics/git/how_to_install_git/)

```
### Install the podman Container package
### tools to install packages in ubuntu box as a root user.
apt-get update -y
#Step 2 – Install Podman
#First, you will need to install some dependencies required to install Podman. You can install them with the following command:

apt-get install curl wget gnupg2 -y
#Next, source your Ubuntu release and add the Podman repository with the following command:

source /etc/os-release
sh -c "echo 'deb http://download.opensuse.org/repositories/devel:/kubic:/libcontainers:/stable/xUbuntu_${VERSION_ID}/ /' > /etc/apt/sources.list.d/devel:kubic:libcontainers:stable.list"
Next, download and add the GPG key with the following command:

wget -nv https://download.opensuse.org/repositories/devel:kubic:libcontainers:stable/xUbuntu_${VERSION_ID}/Release.key -O- | apt-key add -
#Next, update the repository and install Podman with the following command:

apt-get update -qq -y
apt-get -qq --yes install podman
After installing Podman, verify the Podman version with the following command:

podman --version

apt-get install curl wget gnupg2 -y

podman info

If the system is a MACOS do the following to initiate a VM.
#podman machine init -v $HOME:$HOME

```


We will be installing the following tools.

- GenPipes  (https://bitbucket.org/mugqic/genpipes/src/master/)
- Genpipes Container (https://hub.docker.com/r/c3genomics/genpipes)
- Additional bioinformatics tools like TRIMMOMATIC, BOWTIE, TRINITY are available in containers like docker,podman & singularity.


<font color=green>*Note: GenPipes is sponsored by C3G (Canadian center for Computational Genomics)* </font>


<br>

### __2. Workshop files and folders__

Please install the pre-requisite files and softwares before the start of the workshop.


```
## please install Git

git clone https://github.com/dksenthil/denovo_RNAseq_GenPipes_Workshop2022.git


cd denovo_RNAseq_GenPipes_Workshop2022



## Download and uncompress the *tar.gz file.

Google drive link https://drive.google.com/file/d/1nFePjBFOjRuq5lq4_UzgstMYpBTlEegB/view?usp=sharing

## You can use gdrive or wget or browser to download this file. 

tar -xvf rnaseq_CC.tar.qz 


#### File tree


cd $PWD 
> tree -L 3 
.
├── data
│   ├── 01_run_rnaseq_denovoassembly_genepipes.sh
│   ├── 02_rnaseq_genpipesscripts.sh
│   ├── design.rnaseq.txt
│   ├── raw_data
│   │   ├── rnaseq_GM12878_chr19_Rep1_R1.fastq.gz
│   │   ├── rnaseq_GM12878_chr19_Rep1_R2.fastq.gz
│   │   ├── rnaseq_GM12878_chr19_Rep2_R1.fastq.gz
│   │   ├── rnaseq_GM12878_chr19_Rep2_R2.fastq.gz
│   │   ├── rnaseq_H1ESC_chr19_Rep1_R1.fastq.gz
│   │   ├── rnaseq_H1ESC_chr19_Rep1_R2.fastq.gz
│   │   ├── rnaseq_H1ESC_chr19_Rep2_R1.fastq.gz
│   │   └── rnaseq_H1ESC_chr19_Rep2_R2.fastq.gz
│   └── readset.rnaseq.txt
├── data.tar.gz
├── rnaseq_CC
│   └── rnaseq
│       ├── 01_run_rnaseq_denovoassembly_genepipes.sh
│       ├── 02_rnaseq_genpipesscripts.sh
│       ├── 02_rnaseq_genpipesscripts.sh.bk
│       ├── align_and_estimate_abundance
│       ├── blast
│       ├── design.rnaseq.txt
│       ├── differential_expression
│       ├── exploratory
│       ├── filtered_assembly
│       ├── insilico_read_normalization
│       ├── job_output
│       ├── metrics
│       ├── raw_data
│       ├── readset.rnaseq.txt
│       ├── report
│       ├── RnaSeqDeNovoAssembly.20221012T133724.config.trace.ini
│       ├── RnaSeqDeNovoAssembly.20221012T134036.config.trace.ini
│       ├── RnaSeqDeNovoAssembly.20221012T134138.config.trace.ini
│       ├── RnaSeqDeNovoAssembly.20221012T134222.config.trace.ini
│       ├── RnaSeqDeNovoAssembly.20221016T124116.config.trace.ini
│       ├── RnaSeqDeNovoAssembly.config.trace.ini
│       ├── Rplots.pdf
│       ├── submit.log
│       ├── TMHMM_269976
│       ├── trim
│       ├── trinity_out_dir
│       └── trinotate
└── rnaseq_CC.tar.qz

18 directories, 27 files
abacus2:~/C3G_projects/RNASeq_Trinity_Tuxedo_Workshop$


```



<br>


### __3. Resource__
Documentation: https://genpipes.readthedocs.io/en/latest/

GenPipes Container : https://github.com/c3g/genpipes_in_a_container



<br>


### __4. Install Genpipes and dependencies__


```

cd /tmp

### git clone https://bitbucket.org/mugqic/genpipes 

git clone https://bitbucket.org/mugqic/genpipes
Cloning into 'genpipes'...
remote: Enumerating objects: 59845, done.
remote: Counting objects: 100% (3846/3846), done.
remote: Compressing objects: 100% (1008/1008), done.
remote: Total 59845 (delta 3276), reused 3120 (delta 2790), pack-reused 55999
Receiving objects: 100% (59845/59845), 57.72 MiB | 25.31 MiB/s, done.
Resolving deltas: 100% (46666/46666), done.



# move to any folder of interest

mv genpipes /Users/skailasam/Documents/denovo_RNAseq_GenPipes_Workshop2022

```

#### Intialise podman and mount dependency
 
```
### Intialise podman and mount dependency
#>(base) statice:/Users/skailasam/Documents/denovo_RNAseq_GenPipes_Workshop2022
$ podman machine init -v $HOME:$HOME
Error: podman-machine-default: VM already exists

### Python packages 
 podman run  --security-opt label=disable    --rm   --device /dev/fuse --cap-add SYS_ADMIN  -it  -v /tmp:/tmp -v /home/dksenthil:/home/dksenthil    --mount type=volume,source=cvmfs_cache,destination=/cvmfs-cache  c3genomics/genpipes 
 

#>(base) statice:/Users/skailasam/Documents/denovo_RNAseq_GenPipes_Workshop2022
 podman run  --security-opt label=disable    --rm   --device /dev/fuse --cap-add SYS_ADMIN  -it  -v /tmp:/tmp -v /Users/skailasam:/Users/skailasam    --mount type=volume,source=cvmfs_cache,destination=/cvmfs-cache  c3genomics/genpipes

CernVM-FS: running with credentials 999:997
CernVM-FS: loading Fuse module... done
CernVM-FS: mounted cvmfs on /cvmfs/cvmfs-config.computecanada.ca
CernVM-FS: running with credentials 999:997
CernVM-FS: loading Fuse module... done
CernVM-FS: mounted cvmfs on /cvmfs/ref.mugqic
CernVM-FS: running with credentials 999:997
CernVM-FS: loading Fuse module... done
CernVM-FS: mounted cvmfs on /cvmfs/soft.mugqic


```

<br>

### __5. genpipes options and parameters__

```
[root@f5fecc750fec data]# rnaseq_denovo_assembly.py  -h
usage: rnaseq_denovo_assembly.py [-h] [--help] [-c CONFIG [CONFIG ...]] [-s STEPS] [-o OUTPUT_DIR]
                                 [-j {pbs,batch,daemon,slurm}] [-f] [--no-json] [--report] [--clean]
                                 [-l {debug,info,warning,error,critical}] [--sanity-check]
                                 [--force_mem_per_cpu FORCE_MEM_PER_CPU] [--container {wrapper, singularity} <IMAGE PATH>]
                                 [--genpipes_file GENPIPES_FILE] [-t {trinity,seq2fun}] [-d DESIGN] [-r READSETS] [-v]

Version: 4.3.1

For more documentation, visit our website: https://bitbucket.org/mugqic/genpipes/

optional arguments:
  -h                    show this help message and exit
  --help                show detailed description of pipeline and steps
  -c CONFIG [CONFIG ...], --config CONFIG [CONFIG ...]
                        config INI-style list of files; config parameters are overwritten based on files order
  -s STEPS, --steps STEPS
                        step range e.g. '1-5', '3,6,7', '2,4-8'
  -o OUTPUT_DIR, --output-dir OUTPUT_DIR
                        output directory (default: current)
  -j {pbs,batch,daemon,slurm}, --job-scheduler {pbs,batch,daemon,slurm}
                        job scheduler type (default: slurm)
  -f, --force           force creation of jobs even if up to date (default: false)
  --no-json             do not create JSON file per analysed sample to track the analysis status (default: false i.e. JSON
                        file will be created)
  --report              create 'pandoc' command to merge all job markdown report files in the given step range into HTML, if
                        they exist; if --report is set, --job-scheduler, --force, --clean options and job up-to-date status
                        are ignored (default: false)
  --clean               create 'rm' commands for all job removable files in the given step range, if they exist; if --clean
                        is set, --job-scheduler, --force options and job up-to-date status are ignored (default: false)
  -l {debug,info,warning,error,critical}, --log {debug,info,warning,error,critical}
                        log level (default: info)
  --sanity-check        run the pipeline in `sanity check mode` to verify that all the input files needed for the pipeline to
                        run are available on the system (default: false)
  --force_mem_per_cpu FORCE_MEM_PER_CPU
                        Take the mem input in the ini file and force to have a minimum of mem_per_cpu by correcting the
                        number of cpu (default: None)
  --container {wrapper, singularity} <IMAGE PATH>
                        Run inside a container providing a valid singularity image path
  --genpipes_file GENPIPES_FILE, -g GENPIPES_FILE
                        Command file output path. This is the command used to process the data, or said otherwise, this
                        command will "run the Genpipes pipeline". Will be redirected to stdout if the option is not provided.
  -t {trinity,seq2fun}, --type {trinity,seq2fun}
                        Type of pipeline (default trinity)
  -d DESIGN, --design DESIGN
                        design file
  -r READSETS, --readsets READSETS
                        readset file
  -v, --version         show the version information and exit

Steps:
------

----
trinity:
1- picard_sam_to_fastq
2- trimmomatic
3- merge_trimmomatic_stats
4- insilico_read_normalization_readsets
5- insilico_read_normalization_all
6- trinity
7- exonerate_fastasplit
8- blastx_trinity_uniprot
9- blastx_trinity_uniprot_merge
10- transdecoder
11- hmmer
12- rnammer_transcriptome
13- blastp_transdecoder_uniprot
14- signalp
15- tmhmm
16- trinotate
17- align_and_estimate_abundance_prep_reference
18- align_and_estimate_abundance
19- gq_seq_utils_exploratory_analysis_rnaseq_denovo
20- differential_expression
21- filter_annotated_components
22- gq_seq_utils_exploratory_analysis_rnaseq_denovo_filtered
23- differential_expression_filtered
----
seq2fun:
1- picard_sam_to_fastq
2- merge_fastq
3- seq2fun
4- differential_expression_seq2fun
5- pathway_enrichment_seq2fun
```

<br>

### __6. Run Genpipes RNAseq de novo assembly on sample data__

```
[root@f5fecc750fec data]#
rnaseq_denovo_assembly.py --no-json -c $MUGQIC_PIPELINES_HOME/pipelines/rnaseq_denovo_assembly/rnaseq_denovo_assembly.base.ini  -j batch    -d design.rnaseq.txt  -r readset.rnaseq.txt > 02_rnaseq_genpipesscripts.sh


INFO:core.config:Check modules...
INFO:core.config:module check finished

WARNING:core.pipeline:No step provided by the user => launching the entire pipeline

INFO:core.pipeline:Create jobs for step picard_sam_to_fastq...
INFO:core.readset:Parse Illumina readset file readset.rnaseq.txt ...
INFO:core.readset:4 readsets parsed
INFO:core.readset:4 samples parsed

INFO:core.pipeline:Step picard_sam_to_fastq: 0 job created... skipping

INFO:core.pipeline:Create jobs for step trimmomatic...
INFO:core.pipeline:Step trimmomatic: 4 jobs created

INFO:core.pipeline:Create jobs for step merge_trimmomatic_stats...
INFO:core.pipeline:Step merge_trimmomatic_stats: 1 job created

INFO:core.pipeline:Create jobs for step insilico_read_normalization_readsets...
INFO:core.pipeline:Step insilico_read_normalization_readsets: 4 jobs created

INFO:core.pipeline:Create jobs for step insilico_read_normalization_all...
INFO:core.pipeline:Step insilico_read_normalization_all: 2 jobs created

INFO:core.pipeline:Create jobs for step trinity...
INFO:core.pipeline:Step trinity: 2 jobs created

INFO:core.pipeline:Create jobs for step exonerate_fastasplit...
INFO:core.pipeline:Step exonerate_fastasplit: 1 job created

INFO:core.pipeline:Create jobs for step blastx_trinity_uniprot...
INFO:core.pipeline:Step blastx_trinity_uniprot: 20 jobs created

INFO:core.pipeline:Create jobs for step blastx_trinity_uniprot_merge...
INFO:core.pipeline:Step blastx_trinity_uniprot_merge: 2 jobs created

INFO:core.pipeline:Create jobs for step transdecoder...
INFO:core.pipeline:Step transdecoder: 1 job created

INFO:core.pipeline:Create jobs for step hmmer...
INFO:core.pipeline:Step hmmer: 1 job created

FO:core.pipeline:Create jobs for step rnammer_transcriptome...
INFO:core.pipeline:Step rnammer_transcriptome: 1 job created

INFO:core.pipeline:Create jobs for step blastp_transdecoder_uniprot...
INFO:core.pipeline:Step blastp_transdecoder_uniprot: 1 job created

INFO:core.pipeline:Create jobs for step signalp...
INFO:core.pipeline:Step signalp: 1 job created

INFO:core.pipeline:Create jobs for step tmhmm...
INFO:core.pipeline:Step tmhmm: 1 job created

INFO:core.pipeline:Create jobs for step trinotate...
INFO:core.pipeline:Step trinotate: 2 jobs created

INFO:core.pipeline:Create jobs for step align_and_estimate_abundance_prep_reference...
INFO:core.pipeline:Step align_and_estimate_abundance_prep_reference: 1 job created

INFO:core.pipeline:Create jobs for step align_and_estimate_abundance...
INFO:core.pipeline:Step align_and_estimate_abundance: 7 jobs created

INFO:core.pipeline:Create jobs for step gq_seq_utils_exploratory_analysis_rnaseq_denovo...
INFO:core.pipeline:Step gq_seq_utils_exploratory_analysis_rnaseq_denovo: 2 jobs created

INFO:core.pipeline:Create jobs for step differential_expression...
INFO:core.design:Contrast H1ESC_GM12787 (controls: 2, treatments: 2) created
INFO:core.design:1 contrast parsed

INFO:core.pipeline:Step differential_expression: 3 jobs created

INFO:core.pipeline:Create jobs for step filter_annotated_components...
INFO:core.pipeline:Step filter_annotated_components: 2 jobs created

INFO:core.pipeline:Create jobs for step gq_seq_utils_exploratory_analysis_rnaseq_denovo_filtered...
INFO:core.pipeline:Step gq_seq_utils_exploratory_analysis_rnaseq_denovo_filtered: 3 jobs created

INFO:core.pipeline:Create jobs for step differential_expression_filtered...
INFO:core.pipeline:Step differential_expression_filtered: 4 jobs created

INFO:core.pipeline:TOTAL: 66 jobs created

INFO:core.scheduler:
	 To run the script use:
	"  ./<command>.sh"
```

### __5. Run the batch script__

```

bash 02_rnaseq_genpipesscripts.sh

```

<br>

### __6. Output files and format__

```




	
	
