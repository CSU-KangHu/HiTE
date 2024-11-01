FROM continuumio/miniconda3

# Author and maintainer
MAINTAINER Kang Hu <kanghu@csu.edu.cn>
LABEL description="HiTE: A fast and accurate dynamic boundary adjustment approach for full-length Transposable Elements detection and annotation in Genome Assemblies" \
      author="kanghu@csu.edu.cn"

ARG DNAME="HiTE"

RUN apt-get update && apt-get install unzip --yes && apt-get install less --yes && apt-get install curl --yes

# Download HiTE from Github
# RUN git clone https://github.com/CSU-KangHu/HiTE.git

# use local HiTE
COPY HiTE-master.zip /HiTE-master.zip
RUN cd / && unzip HiTE-master.zip && mv HiTE-master /HiTE && rm HiTE-master.zip

# Download HiTE from Zenodo
#RUN curl -LJO https://zenodo.org/records/10574301/files/CSU-KangHu/HiTE-v.3.1.1.zip &&  \
#    unzip HiTE-v.3.1.1.zip && mv CSU-KangHu-HiTE-* /HiTE

RUN conda install mamba -c conda-forge -y
RUN cd /HiTE && chmod +x tools/* bin/NeuralTE/tools/* bin/LTR_FINDER_parallel-master/bin/LTR_FINDER.x86_64-1.0.7/ltr_finder && mamba env create --name ${DNAME} --file=environment.yml && conda clean -a

# Make RUN commands use the new environment
# name need to be the same with the above ${DNAME}
SHELL ["conda", "run", "-n", "HiTE", "/bin/bash", "-c"]

# avoid different perl version conflict
ENV PERL5LIB /
ENV PATH /opt/conda/envs/${DNAME}/bin:$PATH
USER root

WORKDIR /HiTE
RUN cd /HiTE

CMD ["bash"]