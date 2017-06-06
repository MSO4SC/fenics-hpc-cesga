fenics-hpc-cesga
===================

[![Join the chat at https://gitter.im/MSO4SC/fenics-hpc-cesga](https://badges.gitter.im/MSO4SC/fenics-hpc-cesga.svg)](https://gitter.im/MSO4SC/fenics-hpc-cesga?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

Files to compile and execute FEniCS HPC in Finis Terrae II (CESGA). - http://www.fenicsproject.org

GNU and Intel compilers can be used to perform this task. Scripts to use this compilers can be found in **gnu_compiler** and **intel_compiler** folders. This folders are referenced by $COMPILER_DIR variable in the following descriptions.

Install
-------------
 - Login in Finis Terrae II, download and extract fenics-hpc (better in $LUSTRE storage):  
`$ cd $LUSTRE`  
`$ wget http://www.csc.kth.se/~jjan/hpfem2016/fenics-hpc_hpfem.zip`  
`$ unzip fenics-hpc_hpfem.zip`  
 
 - Copy **build_ft2.sh** and **init_ft2.sh** in this repository inside the fenics-hpc_hpfem folder.  
 - Source **build_ft2.sh** to build and install all the FEniCS modules.  
`$ source $COMPILER_DIR/build_ft2.sh`  

 - Build unicorn-minimal  
`$ cd unicorn-minimal`  
`$ make -j 4`  

Execute cube simulation
-------------
For the cube case, the primary aim is to verify that the simulation framework works correctly, and
to explore the adaptiviy. Each adaptive iteration is stored in the folder iter XX with XX the number
of the adaptive iteration.

- Copy **job-adaptive_ft2.script** in this repository to $LUSTRE/fenics-hpc_hpfem/unicorn-minimal/cube_sim01.  
- Execute the demo simulation:  
`$ cd $LUSTRE/fenics-hpc_hpfem`  
`$ source $COMPILER_DIR/init_ft2.sh`  
`$ cd $unicorn-minimal/cube_sim01`  
`$ cp ../demo`  
`$ sbatch job-adaptive_ft2.script`  
