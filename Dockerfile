# syntax=docker/dockerfile:1.4
# Copyright 2023 The University of Hong Kong, Department of Computer Science
# All rights reserved.
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its
#    contributors may be used to endorse or promote products derived from
#    this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

# ClairS Docker image.
#
# Build:
#   docker build -f ./Dockerfile -t hkubal/clairs:v0.5.1 .
#   docker run -it hkubal/clairs:v0.5.1 /opt/bin/run_clairs --help
#
# Build host requires Docker >= 20.10.10. glibc >= 2.34 creates threads with
# clone3, which older Docker seccomp profiles block; glibc does not fall back to
# clone on EPERM, so mamba aborts with "can't start new thread".

FROM condaforge/miniforge3:latest

ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
WORKDIR /opt/bin

# clair3 2.0.3 brings in pytorch-cpu, samtools, whatshap, longphase, parallel,
# pigz, numpy, h5py, hdf5plugin, numexpr, torchmetrics and tqdm, and ships the
# germline models under bin/models. Added here: ClairS runtime deps (scipy,
# scikit-learn, torchinfo) and the toolchain/headers needed to build the
# bundled native components.
#
# libboost-headers rather than boost-cpp: boost-cpp pins icu<=75 and
# zlib 1.2/1.3.1, which cannot co-solve with clair3 2.0.3 (libcurl>=8.22 ->
# libpsl -> icu>=78.3). realign only uses header-only boost.
RUN mamba create -n clairs \
      -c conda-forge \
      -c bioconda \
      python=3.11 \
      clair3=2.0.3 \
      scipy \
      scikit-learn \
      torchinfo \
      gcc \
      gxx \
      make \
      binutils \
      automake \
      libtool \
      libboost-headers \
      curl \
      libcurl \
      zlib \
      bzip2 \
      xz \
      perl \
      -y && \
    mamba clean --all -y

ENV CONDA_PREFIX=/opt/conda/envs/clairs
ENV CONDA_DEFAULT_ENV=clairs
ENV PATH=${CONDA_PREFIX}/bin:/opt/bin:/opt/conda/bin:${PATH}

# pypy is used by several ClairS stages; run_clairs requires pypy >= 3.6.
RUN wget -q https://downloads.python.org/pypy/pypy3.11-v7.3.20-linux64.tar.bz2 && \
    tar -xjf pypy3.11-v7.3.20-linux64.tar.bz2 && \
    rm pypy3.11-v7.3.20-linux64.tar.bz2 && \
    ln -sf /opt/bin/pypy3.11-v7.3.20-linux64/bin/pypy3 ${CONDA_PREFIX}/bin/pypy3 && \
    ln -sf /opt/bin/pypy3.11-v7.3.20-linux64/bin/pypy3 ${CONDA_PREFIX}/bin/pypy && \
    pypy3 -m ensurepip && \
    pypy3 -m pip install --no-cache-dir mpmath==1.2.1

COPY . .

# Native components. All are standalone binaries, independent of the Python
# version: realigner and debruijn_graph from src/realign, and alleleCounter
# (Sanger alleleCount, AGPL) whose setup.sh builds htslib 1.11 and libdeflate
# from source.
RUN cd /opt/bin/src/realign && \
    g++ -std=c++14 -O1 -shared -fPIC -o realigner ssw_cpp.cpp ssw.c realigner.cpp && \
    g++ -std=c++11 -shared -fPIC -o debruijn_graph -O3 debruijn_graph.cpp && \
    cd /opt/bin/src/verdict/allele_counter && \
    chmod +x setup.sh && \
    /bin/bash setup.sh /opt/bin/src/verdict/allele_counter

# Somatic models and CNV reference data. Override the URLs with build args.
ARG CLAIRS_MODELS_URL=https://www.bio8.cs.hku.hk/clairs/models/clairs_models.tar.gz
ARG CLAIRS_CNV_URL=https://www.bio8.cs.hku.hk/clairs/data/reference_files.tar.gz

RUN mkdir -p /opt/models && \
    wget -q "${CLAIRS_MODELS_URL}" -P /opt/models && \
    mkdir -p ${CONDA_PREFIX}/bin/clairs_models && \
    tar -zxf /opt/models/clairs_models.tar.gz -C ${CONDA_PREFIX}/bin/clairs_models && \
    rm /opt/models/clairs_models.tar.gz && \
    mkdir -p /opt/cnv_data ${CONDA_PREFIX}/bin/cnv_data && \
    wget -q "${CLAIRS_CNV_URL}" -P /opt/cnv_data && \
    tar -zxf /opt/cnv_data/reference_files.tar.gz -C ${CONDA_PREFIX}/bin/cnv_data && \
    rm -f /opt/cnv_data/reference_files.tar.gz

# Germline models need no download: the clair3 package ships them under
# bin/models (21 directories, covering every platform run_clairs uses).
# The line below only silences the GNU parallel citation notice.
RUN mkdir -p ~/.parallel && touch ~/.parallel/will-cite

CMD ["/bin/bash"]
