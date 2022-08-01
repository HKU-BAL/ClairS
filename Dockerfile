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
    conda config --add channels anaconda && \
    conda create -n python=3.8 pypy tensorflow>=2.6.0 whatshap samtools -y

ENV PATH /opt/conda/envs/somatic/bin:$PATH
ENV CONDA_DEFAULT_ENV somatic

RUN /bin/bash -c "source activate somatic" && \
    conda install -c conda-forge  -y && \
    conda install -c conda-forge -c bioconda samtools whatshap -y && \
    pip install tensorflow whatshap=1.4 tables==3.6.1 tensorflow-addons==0.16.1 && \
    conda install -c anaconda pigz==2.4 cffi==1.14.4 -y && \
    conda install -c conda-forge parallel=20191122 automake curl xz zlib bzip2 -y && \
    rm -rf /opt/conda/pkgs/* && \
    rm -rf /root/.cache/pip && \

COPY . .
