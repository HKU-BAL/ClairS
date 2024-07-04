# BSD 3-Clause License
#
# Copyright 2023 The University of Hong Kong, Department of Computer Science
# All rights reserved.
#
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

# Example command:
# $ git clone https://github.com/HKU-BAL/ClairS.git
# $ cd ClairS
# $ cd deepvariant
# $ docker build -f ./Dockerfile -t hkubal/clairs:latest .
# $ docker run -it hkubal/clairs:latest /opt/bin/run_clairs --help

FROM ubuntu:16.04

ENV LANG=C.UTF-8 LC_ALL=C.UTF-8 PATH=/opt/bin:/opt/conda/bin:$PATH

# update ubuntu packages
RUN apt-get update --fix-missing && \
    yes|apt-get upgrade && \
    apt-get install -y \
        wget \
        bzip2 \
        make \
        g++ \
        time \
        libboost-graph-dev && \
    rm -rf /bar/lib/apt/lists/*

WORKDIR /opt/bin

# install anaconda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    bash Miniconda3-latest-Linux-x86_64.sh -b -p /opt/conda && \
    rm Miniconda3-latest-Linux-x86_64.sh && \
    conda config --add channels defaults && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge && \
    conda create -n clairs -c pytorch -c conda-forge -c bioconda clair3 pytorch torchinfo tqdm -y && \
    rm -rf /opt/conda/pkgs/* && \
    rm -rf /root/.cache/pip

ENV PATH /opt/conda/envs/clairs/bin:$PATH
ENV CONDA_DEFAULT_ENV clairs

RUN apt install curl zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev -y && \
    /opt/conda/bin/python3 -m pip install scipy scikit-learn && \
    rm -rf /var/lib/apt/lists/*

COPY . .

RUN /bin/bash -c "source activate clairs" && cd /opt/bin/src/realign && \
    g++ -std=c++14 -O1 -shared -fPIC -o realigner ssw_cpp.cpp ssw.c realigner.cpp && \
    g++ -std=c++11 -shared -fPIC -o debruijn_graph -O3 debruijn_graph.cpp && \
    cd /opt/bin/src/verdict/allele_counter && chmod +x setup.sh && /bin/bash setup.sh /opt/bin/src/verdict/allele_counter && \
    wget http://www.bio8.cs.hku.hk/clairs/models/clairs_models.tar.gz	 -P /opt/models && \
    mkdir -p /opt/conda/envs/clairs/bin/clairs_models && \
    tar -zxvf /opt/models/clairs_models.tar.gz -C /opt/conda/envs/clairs/bin/clairs_models && \
    rm /opt/models/clairs_models.tar.gz && \
    mkdir -p /opt/conda/envs/clairs/bin/cnv_data && \
    wget http://www.bio8.cs.hku.hk/clairs/data/reference_files.tar.gz -P /opt/cnv_data && \
    tar -zxvf /opt/cnv_data/reference_files.tar.gz -C /opt/conda/envs/clairs/bin/cnv_data && rm -rf /opt/cnv_data/reference_files.tar.gz && \
    echo 'will cite' | parallel --citation || true \
    echo "source activate clairs" > ~/.bashrc

