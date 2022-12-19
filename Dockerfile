FROM continuumio/miniconda3:latest

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

# install packages
RUN conda config --add channels defaults && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge && \
    conda create -n clair_somatic -c pytorch -c conda-forge -c bioconda clair3 pytorch tqdm torchinfo -y && \
    rm -rf /opt/conda/pkgs/* && \
    rm -rf /root/.cache/pip && \
    echo "source activate clair_somatic" > ~/.bashrc

ENV PATH /opt/conda/envs/clair_somatic/bin:$PATH
ENV CONDA_DEFAULT_ENV clair_somatic

COPY . .
