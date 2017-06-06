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

module purge

module load gcc/5.3.0
module load impi
module load petsc
module load parmetis
module load zlib

PREFIX=$PWD/local
export PKG_CONFIG_PATH=$PREFIX/lib/pkgconfig:$PKG_CONFIG_PATH
export PATH=$PREFIX/bin:$PATH
export PYTHONPATH=$PREFIX/lib64/python2.6/site-packages

wget http://www.csc.kth.se/~jjan/hpfem2016/fenics-hpc_hpfem.zip
unzip fenics-hpc_hpfem.zip
cd fenics-hpc_hpfem

# UFC
cd ufc2-hpc
rm CMakeCache.txt 
cmake -DCMAKE_INSTALL_PREFIX:PATH=$PREFIX
make install
cd ..

# SymPy, Instant, FIAT, UFL, FFC, OrderedDict
for pkg in sympy-0.7.5 instant fiat ufl-1.0.0 ffc-1.0.0 ordereddict-1.1
do
  cd $pkg
  python setup.py install --prefix=$PREFIX
  cd ..
done

source /usr/lib64/xml2Conf.sh

mkdir gts
cd gts
wget pkgs.fedoraproject.org/lookaside/pkgs/gts/gts-snapshot-121130.tar.gz/023ebb6b13b8707534182a3ef0d12908/gts-snapshot-121130.tar.gz
tar xzvf gts-snapshot-121130.tar.gz
cd gts-snapshot-121130
sed -i "/s/Requires: glib-2.0,gthread-2.0,gmodule-2.0/Requires: glib-2.0/g" gts.pc.in
./configure --prefix=$PREFIX --enable-static --disable-shared --with-glib-prefix=$PREFIX --with-pic
make install
make
cd ../../


# DOLFIN-HPC
cd dolfin-hpc
cp /usr/share/aclocal/pkg.m4 m4/
cp /usr/share/aclocal/libxml.m4 m4/
sh regen.sh
CC=gcc CXX=g++ CFLAGS="-O2" CXXFLAGS="-O2" ./configure --prefix=$PREFIX --with-pic --disable-boost-tr1 --with-parmetis --with-petsc=/opt/cesga/petsc/3.7.0/gcc/5.3.0/impi/5.1/ --enable-mpi --enable-mpi-io --disable-progress-bar --disable-xmltest --enable-ufl --with-gts
make -j 8 install
cd ..
cp -av dolfin-hpc/site-packages/dolfin_utils $PREFIX/lib64/python2.6/site-packages
