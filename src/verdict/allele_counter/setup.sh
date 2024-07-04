#!/bin/bash

##########LICENCE##########
# Copyright (c) 2014-2020 Genome Research Ltd.
#
# Author: CASM/Cancer IT <cgphelp@sanger.ac.uk>
#
# This file is part of alleleCount.
#
# alleleCount is free software: you can redistribute it and/or modify it under
# the terms of the GNU Affero General Public License as published by the Free
# Software Foundation; either version 3 of the License, or (at your option) any
# later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more
# details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
##########LICENCE##########

# ALL tool versions used by opt-build.sh
# need to keep in sync with Dockerfile
export VER_HTSLIB="1.11"
export VER_LIBDEFLATE="v1.6"


if [[ ($# -ne 1 && $# -ne 2) ]] ; then
  echo "Please provide an installation path and optionally perl lib paths to allow, e.g."
  echo "  ./setup.sh /opt/myBundle"
  echo "OR all elements versioned:"
  echo "  ./setup.sh /opt/cgpPinel-X.X.X /opt/cgpVcf-X.X.X/lib/perl:/opt/PCAP-core-X.X.X/lib/perl"
  exit 1
fi

INST_PATH=$1

if [[ $# -eq 2 ]] ; then
  CGP_PERLLIBS=$2
fi

# get current directory
INIT_DIR=`pwd`

set -e

# cleanup inst_path
mkdir -p $INST_PATH/bin
cd $INST_PATH
INST_PATH=`pwd`
cd $INIT_DIR

## make sure that build is self contained
#PERLROOT=$INST_PATH/lib/perl5

## allows user to knowingly specify other PERL5LIB areas.
#if [ -z ${CGP_PERLLIBS+x} ]; then
#  export PERL5LIB="$PERLROOT"
#else
#  export PERL5LIB="$PERLROOT:$CGP_PERLLIBS"
#fi

export OPT=$INST_PATH

bash build/opt-build.sh $INST_PATH
bash build/opt-build-local.sh $INST_PATH

echo
echo
echo "Please add the following to beginning of PATH:"
echo "  $INST_PATH/bin"
echo "Please add the following to beginning of LD_LIBRARY_PATH:"
echo "  $INST_PATH/lib"
echo

exit 0
