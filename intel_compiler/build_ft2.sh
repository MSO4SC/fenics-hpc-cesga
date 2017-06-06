##########################
# LOAD MODULES 
##########################

module purge

module load cmake
module load intel/2016
module load impi/5.1
module load petsc/3.7.5-real-opt
module load zlib/1.2.8

##########################
# EXPORT VARIABLES
##########################

export CURR_DIR=$PWD
export FENICS_HPC_SRC_DIR=$CURR_DIR/fenics-hpc_hpfem
#export PREFIX=/opt/cesga/fenics-hpc/0.8.4/intel/2016/impi/5.1.3.223
export PREFIX=$FENICS_HPC_SRC_DIR/build
mkdir -p $PREFIX

export PATH=$PREFIX/bin:$PATH
export LD_LIBRARY_PATH=$PREFIX/lib:$LD_LIBRARY_PATH
export PKG_CONFIG_PATH=$PREFIX/lib/pkgconfig:$PKG_CONFIG_PATH
PYV=`$(which python) -c "import sys;t='{v[0]}.{v[1]}'.format(v=list(sys.version_info[:2]));sys.stdout.write(t)";`
PYP=$PREFIX/lib/python$PYV/site-packages
export PYTHONPATH=$PYP:$PYTHONPATH

export CC=icc 
export CXX=icpc 
export CFLAGS="-O2" 
export CXXFLAGS="-O2" 

##########################
# DOWNLOAD SOURCES
##########################

wget http://www.csc.kth.se/~jjan/hpfem2016/fenics-hpc_hpfem.zip
unzip fenics-hpc_hpfem.zip

##########################
# SymPy, Instant, FIAT, UFL, FFC, OrderedDict
##########################

for pkg in sympy-0.7.5 instant fiat ufl-1.0.0 ffc-1.0.0 ordereddict-1.1
do
  cd $FENICS_HPC_SRC_DIR/$pkg; python setup.py install --prefix=$PREFIX
done

##########################
# UFC
##########################

cd $FENICS_HPC_SRC_DIR/ufc2-hpc
rm CMakeCache.txt 
cmake -DCMAKE_INSTALL_PREFIX:PATH=$PREFIX
make install -j 8

##########################
# GLIB
##########################

cd $FENICS_HPC_SRC_DIR
wget https://ftp.gnome.org/pub/gnome/sources/glib/2.4/glib-2.4.0.tar.gz
tar xzvf glib-2.4.0.tar.gz
cd $FENICS_HPC_SRC_DIR/glib-2.4.0
echo "glib_cv_stack_grows=no" > f.cache
echo "glib_cv_uscore=no" >> f.cache
echo "ac_cv_func_posix_getpwuid_r=yes" >> f.cache
./configure --prefix=$PREFIX --host=x86_64-unknown-linux-gnu --enable-static --disable-shared --cache-file=f.cache --with-pic
make install
cd ..

##########################
# GTS 
##########################
cd $FENICS_HPC_SRC_DIR
wget pkgs.fedoraproject.org/lookaside/pkgs/gts/gts-snapshot-121130.tar.gz/023ebb6b13b8707534182a3ef0d12908/gts-snapshot-121130.tar.gz
tar xzvf gts-snapshot-121130.tar.gz
cd $FENICS_HPC_SRC_DIR/gts-snapshot-121130
sed -i "/s/Requires: glib-2.0,gthread-2.0,gmodule-2.0/Requires: glib-2.0/g" gts.pc.in
./configure --prefix=$PREFIX --enable-static --disable-shared --with-glib-prefix=$PREFIX --with-pic
make install
make


##########################
#DOLFIN-HPC
##########################

cd $FENICS_HPC_SRC_DIR/dolfin-hpc
cp /usr/share/aclocal/pkg.m4 m4/
cp /usr/share/aclocal/libxml.m4 m4/
./regen.sh
./configure --prefix=$PREFIX --with-pic --enable-function-cache --enable-optimize-p1 --disable-boost-tr1 \
                             --with-parmetis --with-petsc=/opt/cesga/petsc/3.7.5/intel/2016/impi/5.1/opt \
                             --enable-openmp --enable-mpi --enable-mpi-io \
                             --disable-progress-bar --disable-xmltest \
                             --enable-ufl --with-gts
make install -j 8

cp -av site-packages/dolfin_utils $PYP

cd $FENICS_HPC_SRC_DIR/unicorn-minimal
make -j 8

cd $CURR_DIR
