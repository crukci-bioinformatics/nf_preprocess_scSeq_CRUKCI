FROM debian:11.7-slim

LABEL authors="Ash Sawle" \
      version="0.1-snapshot" \
      description="Software package/tool dependencies for the CRUK CI preprocess scRNAseq pipeline."

RUN apt-get update && \
    apt-get install -y build-essential curl git procps && \
    apt-get clean

# Install Java
RUN apt-get install -y openjdk-17-jre-headless && \
    apt-get clean

ENV JAVA_HOME=/usr/lib/jvm/jre-17

# Install Cell Ranger
# NOTE THAT THE URL IS A TEMPORARY LINK AND WILL EXPIRE

RUN curl -o cellranger-8.0.1.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-8.0.1.tar.gz?Expires=1719975091&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA&Signature=LXri15y4wRzFSEHkL~nsAS8P04iVoOaEiml4uBid9GtPQY53iy56vsc~TepB~oTqKTUxoX4f0dGJtmgMnHWt7PhhKxzSUZHpN0boWuLHjgfTAxggNzmfv3fY9YZ1AH7LgCawoUJBJd9jYG32rEE1FdkbT32OADDcvgvAxWYQIGsxR9Eg-irlZcDhQSkW4xhC2Up-7zlZaB5XCnDgIQBSOTzfdrikh5iY3g34fasYtlMBaRjYNyAwrF6bUeOOe2rYZdSpEltoduegwyzl4iq0Tww9yi2cFT6ClktKPzONZLzfZvbKbGG3pp78t1WZ5K-AtjulMOzdPhfZSF3vL~V7hA__" && \
    mkdir -p /opt && \
    tar -C /opt -xzvf cellranger-8.0.1.tar.gz && \
    mv /opt/cellranger-8.0.1 /opt/cellranger && \
    rm -f cellranger-8.0.1.tar.gz

ENV PATH=/opt/cellranger/bin:$PATH

# Install Miniconda
RUN curl -L -O https://github.com/conda-forge/miniforge/releases/download/24.3.0-0/miniforge3-linux-x86_64.sh && \
    mkdir -p /opt && \
    sh miniforge3-linux-x86_64.sh -b -p /opt/miniforge3 && \
    rm miniforge3-linux-x86_64.sh

# Install R using Miniconda

COPY conda_R.yaml .

RUN /opt/miniforge3/bin/mamba env create -f conda_R.yaml && /opt/miniforge3/bin/mamba clean -a

ENV PATH=/opt/miniforge3/envs/R/bin:$PATH

