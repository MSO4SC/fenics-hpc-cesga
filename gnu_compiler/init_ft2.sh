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

source /usr/lib64/xml2Conf.sh
