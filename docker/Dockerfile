FROM ubuntu:14.04
MAINTAINER murray.cadzow@otago.ac.nz
RUN apt-get update && apt-get install -y \
    python-setuptools \
    python-numpy \
    python-scipy \
    git \
    wget

RUN apt-get install -y --no-install-recommends \
    r-base r-base-dev r-recommended littler
RUN ln -s /usr/share/doc/littler/examples/install.r \
    /usr/local/bin/install.r

RUN git clone -b selectionTools1.1-dev https://github.com/smilefreak/selectionTools.git  selectionTools1.1-dev\
    && cd selectionTools1.1-dev \
    && ./install.sh
RUN /usr/local/bin/selection_pipeline -h


