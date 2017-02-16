#!/bin/bash -l
# The -l above is required to get the full environment with modules

# Copyright 2017 MSO4SC - javier.carnero@atos.net
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#   http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# The name of the script is fenics-hpc-orig
#SBATCH -J fenics-cube

# Only 1 hour wall-clock time will be given to this job
#SBATCH -t 01:00:00

# Number of nodes
#SBATCH -N 4
# Number of MPI processes per node (the following is actually the default)
#SBATCH --ntasks-per-node=24
# Number of MPI processes.
#SBATCH -n 96
# Partition
#SBATCH -p thinnodes

#SBATCH -e error_file.e
#SBATCH -o output_file.o

export ATP_ENABLED=1

set -e

# mesh0.bin is the starting mesh
cp mesh0.bin mesh.bin

M=24
MMAX=96

for i in $( seq -w 00 40 )
do
    # Run one adaptive iteration
    mkdir -p iter_${i}
#    ( aprun -n $M -N 32 -ss -j1 ./demo > log1 2> log2 < /dev/null )
    ( srun -n $M ./demo > log1 2> log2 < /dev/null )

    # Postprocess the solution into ParaView files
    srun -n 1 ./dolfin_post -m mesh_out.bin -t vtk -n 200 -s velocity -f 10  1> log.pp1 2> log.ppe1 &
    srun -n 1 ./dolfin_post -m mesh_out.bin -t vtk -n 1000 -s dvelocity -f 10  1> log.pp2 2> log.ppe2 &
    srun -n 1 ./dolfin_post -m mesh_out.bin -t vtk -n 1000 -s pressure -f 10 1> log.pp3 2> log.ppe3 &
    srun -n 1 ./dolfin_post -m mesh_out.bin -t vtk -n 1000 -s dpressure -f 10 1> log.pp4 2> log.ppe4 &
    wait

    # Prepare for the next adaptive iteration
    mv *.bin log1 log2 *.vtu *.pvd iter_${i}
    cp iter_${i}/rmesh.bin iter_${i}/mesh0.bin iter_${i}/log1 .
    cp rmesh.bin mesh.bin
    
    # Compute new core count
    M=$(( $( tail -n 10 log1|grep "vertices after"|cut -d " " -f 3 ) / 250 ))
    M=$(( $M < 24 ? 24 : $M ))
    M=$(( $M > $MMAX ? $MMAX : $M ))
done
