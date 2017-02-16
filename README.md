fenics-hpc-cesga
===================
Files to compile and execute FEniCS HPC in Finis Terrae II (CESGA). - http://www.fenicsproject.org

Install
-------------
 - Login in Finis Terrae II, download and extract fenics-hpc (better in $LUSTRE storage):
 `$ cd $LUSTRE`
 `$ wget http://www.csc.kth.se/~jjan/hpfem2016/fenics-hpc_hpfem.zip`
 `$ unzip fenics-hpc_hpfem.zip`
 
 - Copy **build_ft2.sh** and **init_ft2.sh** in this repository inside the fenics-hpc_hpfem folder
 - Source **build_ft2.sh** to build and install all the FEniCS modules.
 `$ source build_ft2.sh`
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
`$ source init_ft2.sh`
`$ cd $unicorn-minimal/cube_sim01`
`$ cp ../demo`
`$ sbatch job-adaptive_ft2.script`