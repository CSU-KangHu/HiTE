FROM continuumio/miniconda3

# Author and maintainer
MAINTAINER Kang Hu <kanghu@csu.edu.cn>
LABEL description="HiTE: An accurate dynamic boundary adjustment approach for full-length Transposable Elements detection and annotation in Genome Assemblies" \
      author="kanghu@csu.edu.cn"

ARG DNAME="HiTE"

RUN apt-get update && apt-get install unzip --yes

# Command 'RUN' during docker build
#download RepeatMasker libraries
RUN git clone https://github.com/CSU-KangHu/TE_annotation.git && \
    cd TE_annotation && unzip RepeatMasker_Lib.zip

# Download HiTE and create the environment
#COPY environment.yml /
RUN git clone https://github.com/CSU-KangHu/HiTE.git
RUN cd /HiTE && conda env create --name ${DNAME} --file=environment.yml && conda clean -a

#update RepeatMasker libraries
RUN mv /TE_annotation/RepeatMasker_Lib/* /opt/conda/envs/HiTE/share/RepeatMasker/Libraries/


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