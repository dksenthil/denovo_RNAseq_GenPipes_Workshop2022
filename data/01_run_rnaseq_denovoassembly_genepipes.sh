#mv 02_rnaseq_genpipesscripts.sh 02_rnaseq_genpipesscripts.sh.bk

#module purge
#module load mugqic/python/2.7.14
#module load mugqic/genpipes/3.6.1

#local machine 
rnaseq_denovo_assembly.py --no-json  -c $MUGQIC_PIPELINES_HOME/pipelines/rnaseq_denovo_assembly/rnaseq_denovo_assembly.base.ini  -j batch    -d design.rnaseq.txt  -r readset.rnaseq.txt     > 02_rnaseq_genpipesscripts.sh  


