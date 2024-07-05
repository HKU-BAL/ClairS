#! /bin/bash

set -xe

if [[ -z "${TMPDIR}" ]]; then
  TMPDIR=/tmp
fi

set -u

if [ "$#" -lt "1" ] ; then
  echo "Please provide an installation path such as /opt/ICGC"
  exit 1
fi

# get path to this script
SCRIPT_PATH=`dirname $0`;
SCRIPT_PATH=`(cd $SCRIPT_PATH && pwd)`

# get the location to install to
INST_PATH=$1
mkdir -p $1
INST_PATH=`(cd $1 && pwd)`
echo $INST_PATH

# get current directory
INIT_DIR=`pwd`

CPU=`grep -c ^processor /proc/cpuinfo`
if [ $? -eq 0 ]; then
  if [ "$CPU" -gt "6" ]; then
    CPU=6
  fi
else
  CPU=1
fi
echo "Max compilation CPUs set to $CPU"

SETUP_DIR=$INIT_DIR/install_tmp
mkdir -p $SETUP_DIR/distro # don't delete the actual distro directory until the very end
mkdir -p $INST_PATH/bin
cd $SETUP_DIR

# make sure tools installed can see the install loc of libraries
set +u
#export LD_LIBRARY_PATH=`echo $INST_PATH/lib:$LD_LIBRARY_PATH | perl -pe 's/:\$//;'`
#export PATH=`echo $INST_PATH/bin:$PATH | perl -pe 's/:\$//;'`
export LD_LIBRARY_PATH="$INST_PATH/lib:$LD_LIBRARY_PATH"
export PATH="$INST_PATH/bin:$PATH"

#export MANPATH=`echo $INST_PATH/man:$INST_PATH/share/man:$MANPATH | perl -pe 's/:\$//;'`
#export PERL5LIB=`echo $INST_PATH/lib/perl5:$PERL5LIB | perl -pe 's/:\$//;'`
#CPANM=`which cpanm`
#echo  "Installing Perl prerequisites ..."
#$CPANM --no-interactive --notest --mirror http://cpan.metacpan.org -l $INST_PATH/ --installdeps $INIT_DIR/perl/. < /dev/null

set -u
### alleleCount
echo "Building alleleCounter ..."
if [ ! -e $SETUP_DIR/alleleCount.success ]; then
  #build the C part
  cd $INIT_DIR
  mkdir -p $INIT_DIR/c/bin
#  make -C c clean
  export prefix=$INST_PATH
  make -C c -j$CPU
  cp $INIT_DIR/c/bin/alleleCounter $INST_PATH/bin/.
  #build the perl part
#  cd $INIT_DIR/perl
#  perl Makefile.PL INSTALL_BASE=$INST_PATH
#  make
#  make test
#  make install
  touch $SETUP_DIR/alleleCounter.success
fi
