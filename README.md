## __denovo_RNAseq_GenPipes_Workshop2022__

___

 Demo data a, slides and protocol to run Trinity based denovo RNA assembly pipeline in GenPipes. 


___

Author: Senthilkumar Kailasam

Date: Oct 2022


Workshop files
```
git clone https://github.com/dksenthil/denovo_RNAseq_GenPipes_Workshop2022.git

```

We  will go through the steps involved in installing Genpipes in Ubuntu 22.04 LTS


Tools required:

- podman
- bash
- git
- GenPipes
- Additional tools like TRIMMOTIC, BOWTIE2, TRINITY are available in containers like docker,podman,singularity.



<font color=green>*Note: GenPipes is psonsored by C3G (Canadian center for Computational Genomics)* </font>


#### __1. Install all  prerequisites/dependency__

*  Assuming its a fresh installation of Ubuntu  LTS

Documentation: https://genpipes.readthedocs.io/en/latest/


```r bash

sudo apt install git-all

git clone https://bitbucket.org/mugqic/genpipes $HOME/some/dir/genpipes

#Prerequisites: GIT 
```


```python
### commands to install packages in ubuntu box as a root user.
podman machine init -v $HOME:$HOME


### Python packages 
 podman run  --security-opt label=disable    --rm   --device /dev/fuse --cap-add SYS_ADMIN  -it  -v /tmp:/tmp -v /Users/skailasam:/Users/skailasam    --mount type=volume,source=cvmfs_cache,destination=/cvmfs-cache  c3genomics/genpipes 
```

#### __2. GNU GCC 6 Installation__

#### __4. CUDA 9.0 Installation__

* Assuming that the NVIDIA DRIVER installation went well! If not see that instruction at [Step 6](#NVIDIA-DRIVER)
* Download cuda_9.0.176_384.81_linux.run  from NVIDIA website.




#### __5. Compiling AMBER 16 and AMBERTOOLS 17__

* Download AMBER with license purchase at [ambermd.org](http://ambermd.org/GetAmber.php)
* Download AmberTools17 from [ambermd.org](http://ambermd.org/AmberTools17-get.html)




#### __6. NVIDIA DRIVER INSTALLATION__

