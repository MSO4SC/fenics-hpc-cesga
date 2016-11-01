Unicorn
-------

Description
===========

Unicorn is the adaptive continuum mechanics solver of the FEniCS-HPC
branch in the FEniCS project (http://fenicsproject.org), focusing
mainly on turbulent flow and fluid-structure interaction (FSI).

Unicorn stands for Unified Continuum modeling, where we pose the
fundamental balance laws for mass, momentum and energy with a general
Cauchy stress, and allow different constitutive laws chosen by a phase
marker, for example for a Newtonian fluid and a Neo-Hookean
solid.. The total model is then solved using the G2 General
Galerkin/Direct FEM methodology which is based on stabilized FEM and
goal-oriented adaptive error control.

Unicorn is part of the software distribution of FEniCS-HPC
(https://bitbucket.org/fenics-hpc) with automatic installation scripts
for the Beskow supercomputer at KTH. This branch has been developed
mainly by the Computational Technology Laboratory at KTH
(http://ctl.csc.kth.se).

Installation
============

1. Login to Beskow
2. cd to fenics-hpc distribution directory
3. source build_beskow.sh (use: source init_beskow.sh for just setting up the environment)

Testing
=======

1. cd unicorn
2. make -j 4 (parallel build with 4 processes)
3. cd cube_sim01 (directory for example simulation)
4. cp ../demo .  (copy the just-compiled program)
5. sbatch job-adaptive.script (run the adaptive simulation)
6. sh output.sh (to extract and print the key output valies of the simulation)

To visualize, copy the *.pvd, *.pvtu and *.vtu files to your local computer and open the *.pvd files with ParaView.

Contact (for this distribution)
=======

Johan Jansson (jjan@kth.se)
Niclas Jansson (njansson@kth.se)

References
==========

[1] Johan Hoffman, Johan Jansson, Niclas Jansson, FEniCS-HPC:
Automated predictive high-performance finite element computing with
applications in aerodynamics, Proceedings of the 11th International
Conference on Parallel Processing and Applied Mathematics, PPAM
2015. Lecture Notes in Computer Science, 2015
