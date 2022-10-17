## __denovo_RNAseq_GenPipes_Workshop2022__

___

 Demo data a, slides and protocol to run Trinity based denovo RNA assembly pipeline in GenPipes. 


___

Author: Senthilkumar Kailasam

Date: Oct 2022


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

*  Ubuntu 18.04 LTS comes with GCC 7 and AMBER16 requires GCC 6. 


```python
### GCC6 installation
aptitude install gcc-6 g++-6 gfortran-6

#### create Symbolic Links in /usr/bin
rm -f gcc gcc-ar gcc-nm gcc-ranlin g++ gfortran
ln -s gcc-6 gcc
ln -s gcc-ar-6 gcc-ar
ln -s gcc-nm-6 gcc-nm
ln -s gcc-ranlib-6 gcc-ranlib
ln -s g++-6 g++
ln -s /usr/bin/gfortran-6 /usr/bin/gfortran 
```

#### __3. NETCDF installation__

* Download parallel-netcdf-1.9.0 




```python
### Install netcdf

## go to netcdf folder followed by make and make install
cd parallel-netcdf-1.9.0/ 
./configure

make 

make install 
```

#### __4. CUDA 9.0 Installation__

* Assuming that the NVIDIA DRIVER installation went well! If not see that instruction at [Step 6](#NVIDIA-DRIVER)
* Download cuda_9.0.176_384.81_linux.run  from NVIDIA website.


```python
### On the terminal , install cuda
wget https://developer.nvidia.com/compute/cuda/9.0/Prod/local_installers/cuda_9.0.176_384.81_linux-run

./cuda_9.0.176_384.81_linux.run 

### add lines to ld.so.conf

echo "include /usr/local/cuda-9.0/lib64" >> /etc/ld.so.conf 

### run ldconfig

ldconfig 

### Add PATH and environment updates to ~/.bashrc file and source

export PATH="/usr/software/anaconda2/bin:$PATH"
test -f /usr/software/amber16/amber.sh && source /usr/software/amber16/amber.sh

export PATH=/usr/local/cuda-9.0/bin${PATH:+:${PATH}}export 
LD_LIBRARY_PATH=/usr/local/cuda/lib64${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}} 

source ~/.bashrc 
```

#### __5. Compiling AMBER 16 and AMBERTOOLS 17__

* Download AMBER with license purchase at [ambermd.org](http://ambermd.org/GetAmber.php)
* Download AmberTools17 from [ambermd.org](http://ambermd.org/AmberTools17-get.html)


```python
### Steps on the terminal

### Steps to compile the serial verison assuming my files are in tmp folder.

tar xvvf /tmp/Amber16.tar.bz2 

tar xvvf /tmp/AmberTools17.tar.bz2 

cd /usr/software/amber16 

./configure gnu 

make install 

### Steps to compile the parallel verison with NETCDF capabilities

./configure -mpi --with-pnetcdf /usr/software/parallel-netcdf gnu 

make install 

#### Steps to compile serial cuda enabled amber 

./configure -cuda gnu 
make install 
    
#### Steps to compile parallel+netcdf+cuda enabled amber 
./configure -mpi -cuda --with-pnetcdf  /usr/soft/parallel-netcdf gnu 

make install 



```



#### __6. NVIDIA DRIVER INSTALLATION__


I followed the instruction given at 
https://www.mvps.net/docs/install-nvidia-drivers-ubuntu-18-04-lts-bionic-beaver-linux/
<a id='NVIDIA-DRIVER'>
</a>


```python
### Clean the system of other Nvidia drivers

aptitude purge nvidia*

### check the latest driver version for our Nvidia GPU

using link https://www.nvidia.com/object/unix.html
    
### Add the Nvidia graphic card PPA

add-apt-repository ppa:graphics-drivers
    
### Prepare the system for the installation

aptitude update

### Install the Nvidia GPU driver

aptitude install nvidia-driver-390

reboot

#### Verify Nvidia Driver installation

lsmod | grep nvidia
or
nvidia-smi

```
