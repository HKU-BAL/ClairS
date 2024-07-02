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
export LD_LIBRARY_PATH=`echo $INST_PATH/lib:$LD_LIBRARY_PATH | perl -pe 's/:\$//;'`
export PATH=`echo $INST_PATH/bin:$PATH | perl -pe 's/:\$//;'`

export LD_LIBRARY_PATH="$INST_PATH/lib:$LD_LIBRARY_PATH"
export PATH="$INST_PATH/bin:$PATH"

#export MANPATH=`echo $INST_PATH/man:$INST_PATH/share/man:$MANPATH | perl -pe 's/:\$//;'`
#export PERL5LIB=`echo $INST_PATH/lib/perl5:$PERL5LIB | perl -pe 's/:\$//;'`
set -u

### INSTALL CPANMINUS
#curl -sSL https://cpanmin.us/ > $SETUP_DIR/cpanm
#perl $SETUP_DIR/cpanm --no-wget --no-interactive --notest --mirror http://cpan.metacpan.org -l $INST_PATH App::cpanminus
#rm -f $SETUP_DIR/cpanm
#
#echo "Installing Perl base deps ..."
#if [ ! -e $SETUP_DIR/basePerlDeps.success ]; then
#  perlmods=( "ExtUtils::CBuilder" "Module::Build~0.42" "Const::Fast" "File::Which" "LWP::UserAgent")
#  for i in "${perlmods[@]}" ; do
#    cpanm --no-interactive --notest --mirror http://cpan.metacpan.org -l $INST_PATH $i
#  done
#  touch $SETUP_DIR/basePerlDeps.success
#fi

## libdeflate
if [ ! -e $SETUP_DIR/libdeflate.success ]; then
  rm -rf tmp_deflate
  mkdir -p tmp_deflate
  curl -sSL --retry 10 https://github.com/ebiggers/libdeflate/archive/${VER_LIBDEFLATE}.tar.gz > distro.tar.gz
  tar --strip-components 1 -C tmp_deflate -zxf distro.tar.gz
  cd tmp_deflate
  PREFIX=$INST_PATH make -j$CPU CFLAGS="-fPIC -O3" install
  cd ../
  rm -rf distro.* tmp_deflate
  touch $SETUP_DIR/libdeflate.success
fi

SOURCE_HTSLIB="https://github.com/samtools/htslib/releases/download/${VER_HTSLIB}/htslib-${VER_HTSLIB}.tar.bz2"

cd $SETUP_DIR

echo "Downloading htslib ..."
if [ ! -e $SETUP_DIR/htslibGet.success ]; then
  cd $SETUP_DIR
  wget $SOURCE_HTSLIB
  touch $SETUP_DIR/htslibGet.success
fi

echo "Building htslib ..."
if [ ! -e $SETUP_DIR/htslib.success ]; then
  mkdir -p htslib
  tar --strip-components 1 -C htslib -jxf htslib-${VER_HTSLIB}.tar.bz2
  cd htslib
  ./configure --enable-plugins --enable-libcurl --with-libdeflate --prefix=$INST_PATH \
  CPPFLAGS="-I$INST_PATH/include" \
  LDFLAGS="-L${INST_PATH}/lib -Wl,-R${INST_PATH}/lib"
  make -j$CPU
  make install
  cd $SETUP_DIR
  touch $SETUP_DIR/htslib.success
fi
