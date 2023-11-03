FROM continuumio/miniconda3

# Author and maintainer
MAINTAINER Kang Hu <kanghu@csu.edu.cn>
LABEL description="HiTE: A fast and accurate dynamic boundary adjustment approach for full-length Transposable Elements detection and annotation in Genome Assemblies" \
      author="kanghu@csu.edu.cn"

ARG DNAME="HiTE"

RUN apt-get update && apt-get install unzip --yes && apt-get install less --yes && apt-get install curl --yes

# Command 'RUN' during docker build
#download RepeatMasker libraries from Github
#RUN git clone https://github.com/CSU-KangHu/TE_annotation.git && \
#    cd TE_annotation && unzip RepeatMasker_Lib.zip \

#download RepeatMasker libraries from Zenodo
RUN curl -LJO https://zenodo.org/records/10068148/files/CSU-KangHu/TE_annotation-v3.0.zip?download=1 &&  \
    unzip TE_annotation-v3.0.zip && cd CSU-KangHu-TE_annotation-* && unzip RepeatMasker_Lib.zip

# Download HiTE and create the environment
# COPY environment.yml /
# Download HiTE from Github
# RUN git clone https://github.com/CSU-KangHu/HiTE.git
# RUN cd /HiTE && conda env create --name ${DNAME} --file=environment.yml && conda clean -a
# Download HiTE from Zenodo
RUN curl -LJO https://zenodo.org/records/10068550/files/CSU-KangHu/HiTE-v.3.0.1.zip?download=1 &&  \
    unzip HiTE-v.3.0.1.zip && cd CSU-KangHu-HiTE-* && chmod +x tools/*  \
    && conda env create --name ${DNAME} --file=environment.yml && conda clean -a

# update RepeatMasker libraries from Github
# RUN mv /TE_annotation/RepeatMasker_Lib/* /opt/conda/envs/HiTE/share/RepeatMasker/Libraries/
# update RepeatMasker libraries from Zenodo
RUN mv /CSU-KangHu-TE_annotation-*/RepeatMasker_Lib/* /opt/conda/envs/HiTE/share/RepeatMasker/Libraries/

# Make RUN commands use the new environment
# name need to be the same with the above ${DNAME}
SHELL ["conda", "run", "-n", "HiTE", "/bin/bash", "-c"]

# avoid different perl version conflict
ENV PERL5LIB /
ENV PATH /opt/conda/envs/${DNAME}/bin:$PATH
USER root

# HiTE download from Github
# WORKDIR /HiTE
# RUN cd /HiTE

# HiTE download from Zenodo
WORKDIR /CSU-KangHu-HiTE-*
RUN cd /CSU-KangHu-HiTE-*

CMD ["bash"]