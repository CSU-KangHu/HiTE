FROM continuumio/miniconda3

# Author and maintainer
MAINTAINER Kang Hu <kanghu@csu.edu.cn>
LABEL description="HiTE: A progressive method for accurate detection of intact Transposable Elements" \
      author="kanghu@csu.edu.cn"

ARG DNAME="HiTE"

# Command 'RUN' during docker build
# Download HiTE and create the environment
COPY environment.yml /
RUN git clone https://github.com/CSU-KangHu/HiTE.git
RUN conda env create --name ${DNAME} --file=environment.yml && conda clean -a

# Make RUN commands use the new environment
# name need to be the same with the above ${DNAME}
SHELL ["conda", "run", "-n", "HiTE", "/bin/bash", "-c"]

ENV PATH /opt/conda/envs/${DNAME}/bin:$PATH
USER root

WORKDIR /HiTE/ReferenceMode

RUN cd /HiTE/ReferenceMode

CMD ["bash"]