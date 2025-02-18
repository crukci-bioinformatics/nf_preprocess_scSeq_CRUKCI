FROM ubuntu:24.04

LABEL authors="Ash Sawle" \
      version="1.0.0" \
      description="Software package/tool dependencies for the CRUK CI preprocess scRNAseq pipeline."

RUN apt-get update && \
    apt-get install -y build-essential curl git procps && \
    apt-get clean

# Install Java
RUN apt-get install -y openjdk-21-jre-headless && \
    apt-get clean

ENV JAVA_HOME=/usr/lib/jvm/jre-21

# Install Cell Ranger
# NOTE THAT THE URL IS A TEMPORARY LINK AND WILL EXPIRE

RUN curl -o cellranger-9.0.1.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-9.0.1.tar.gz?Expires=1739844570&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA&Signature=DQZJQNqFd1FS1mxa7DTEf3k4sGQlAMU2tnBDmX2v-Y2xOMcGrGzQ0eOJ6p1CArAIBugkuOJhU30tdhXgqxQ8xdL0w5xipirTkOcnRyGDw~i~s3Aw6Cmrab3lPSqbCE5sfW~ZhAAsJE2z9nDSxwoc4~YALVY~If6YQflB3bPYwmJ5a2BcVxxXPgfr4ZBTJV~Q3oCHyJ8ZXgwXOZiKXg2rLgLdddZDNd~UWfcKeg8gbhpD~OY2ZIaWM6vmGKwu9c3uIlg0Chia9VApeTf78uyuZaeI5m8eUQhKg0I8hb0ZD4f-2sMhnjjIUgjsLDNWdfbDa6mjI8eN4ngxseOKiAmngg__" && \
    mkdir -p /opt && \
    tar -C /opt -xzvf cellranger-9.0.1.tar.gz && \
    mv /opt/cellranger-9.0.1 /opt/cellranger && \
    rm -f cellranger-9.0.1.tar.gz

ENV PATH=/opt/cellranger/bin:$PATH

# Install Miniconda
RUN curl -L -O https://github.com/conda-forge/miniforge/releases/download/25.1.1-0/miniforge3-linux-x86_64.sh && \
    mkdir -p /opt && \
    sh miniforge3-linux-x86_64.sh -b -p /opt/miniforge3 && \
    rm miniforge3-linux-x86_64.sh

# Install R using Miniconda

COPY conda_R.yaml .

RUN /opt/miniforge3/bin/mamba env create -f conda_R.yaml && /opt/miniforge3/bin/mamba clean -a

ENV PATH=/opt/miniforge3/envs/R/bin:$PATH


