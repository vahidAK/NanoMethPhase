FROM ubuntu:focal
# NanoMethPhase & SNVoter
ARG DEBIAN_FRONTEND=noninteractive
RUN apt update && apt install -y wget gcc git libz-dev build-essential \
    dirmngr apt-transport-https ca-certificates software-properties-common \
    python3 python3-pip make python-dev libhdf5-dev libcudart10.1
RUN pip3 install --upgrade pip setuptools wheel
RUN pip3 install --use-feature=2020-resolver nanomethphase
RUN apt install -y --no-install-recommends r-base r-cran-biocmanager
RUN Rscript -e "BiocManager::install('DSS')"
# NanoPolish
RUN git clone --recursive https://github.com/jts/nanopolish.git
RUN pip install -r nanopolish/scripts/requirements.txt
RUN cd nanopolish && make && cd ..
# Clair
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
RUN bash Miniconda3-latest-Linux-x86_64.sh -b && \
    rm Miniconda3-latest-Linux-x86_64.sh
ENV PATH="/root/miniconda3/bin:${PATH}"
ARG PATH="/root/miniconda3/bin:${PATH}"
RUN conda create -n clair-env -c conda-forge -c bioconda -y clair
RUN /root/miniconda3/envs/clair-env/bin/pypy3 -m ensurepip && \
    /root/miniconda3/envs/clair-env/bin/pypy3 -m pip install \
    --no-cache-dir intervaltree
RUN mkdir ont
RUN cd ont && wget http://www.bio8.cs.hku.hk/clair_models/ont/122HD34.tar && \
    tar -xf 122HD34.tar && cd ../
# WhatsHap
RUN pip3 install whatshap
# NanoPolish & Clair in path
RUN ln -rs /nanopolish/nanopolish /bin/nanopolish
RUN ln -rs /root/miniconda3/envs/clair-env/bin/clair.py /bin/clair

