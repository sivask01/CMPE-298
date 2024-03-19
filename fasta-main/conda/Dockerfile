# Copyright (c) David C. Anastasiu & Gheorghi Guzun
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

FROM ubuntu

RUN apt-get update && apt-get install -y wget git

RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
        bash Miniconda3-latest-Linux-x86_64.sh -b -p ~/miniconda3
ENV PATH="/root/miniconda3/condabin:${PATH}"

RUN conda install conda-build

COPY ./ fasta
WORKDIR /fasta/conda

RUN conda build fasta --no-anaconda-upload -c pytorch
