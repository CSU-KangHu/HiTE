FROM continuumio/miniconda3

# Author and maintainer
MAINTAINER Kang Hu <kanghu@csu.edu.cn>
LABEL description="HiTE: A progressive method for accurate detection of intact Transposable Elements" \
      author="kanghu@csu.edu.cn"

ARG DNAME="HiTE"

RUN apt-get update && apt-get install unzip --yes

## update aliyun source
#RUN echo 'deb https://mirrors.aliyun.com/debian/ bullseye main non-free contrib' > /etc/apt/sources.list && \
#    echo 'deb-src https://mirrors.aliyun.com/debian/ bullseye main non-free contrib' >> /etc/apt/sources.list && \
#    echo 'deb https://mirrors.aliyun.com/debian-security/ bullseye-security main' >> /etc/apt/sources.list && \
#    echo 'deb-src https://mirrors.aliyun.com/debian-security/ bullseye-security main' >> /etc/apt/sources.list && \
#    echo 'deb https://mirrors.aliyun.com/debian/ bullseye-updates main non-free contrib' >> /etc/apt/sources.list && \
#    echo 'deb-src https://mirrors.aliyun.com/debian/ bullseye-updates main non-free contrib' >> /etc/apt/sources.list && \
#    echo 'deb https://mirrors.aliyun.com/debian/ bullseye-backports main non-free contrib' >> /etc/apt/sources.list && \
#    echo 'deb-src https://mirrors.aliyun.com/debian/ bullseye-backports main non-free contrib' >> /etc/apt/sources.list && \
#    apt-get update
#
##solve git clone SSL error
#RUN apt-get install build-essential fakeroot dpkg-dev -y && \
#    apt-get build-dep git -y && \
#    apt-get install libcurl4-openssl-dev -y && \
#    cd ~ && \
#    mkdir source-git && \
#    cd source-git/ && \
#    apt-get source git && \
#    cd git-2.*.*/ && \
#    sed -i -- 's/libcurl4-gnutls-dev/libcurl4-openssl-dev/' ./debian/control && \
#    sed -i -- '/TEST\s*=\s*test/d' ./debian/rules && \
#    dpkg-buildpackage -rfakeroot -b -uc -us && \
#    dpkg -i ../git_*ubuntu*.deb



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

ENV PERL5LIB /
ENV PATH /opt/conda/envs/${DNAME}/bin:$PATH
USER root

WORKDIR /HiTE/ReferenceMode

RUN cd /HiTE/ReferenceMode

CMD ["bash"]