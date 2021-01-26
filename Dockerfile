FROM ubuntu:focal
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8 PATH=/opt/clair/bin:/opt/conda/bin:$PATH
# apt update & dependencies install
ARG DEBIAN_FRONTEND=noninteractive
RUN echo "deb http://security.ubuntu.com/ubuntu xenial-security main" \
    >> /etc/apt/sources.list \
    && apt update && apt upgrade -y \
    && apt install -y --no-install-recommends \
    wget gcc git libz-dev build-essential dirmngr apt-transport-https \
    ca-certificates software-properties-common python3 python3-pip make \
    python-dev libhdf5-dev libcudart10.1 libssl1.0.0 \
    && apt-get clean
# Clair dir
WORKDIR /opt/clair
COPY entrypoint.sh entrypoint.sh
# NanoMethPhase & SNVoter
RUN pip3 install nanomethphase \
# WhatsHap
    && pip3 install whatshap \
# NanoPolish
    && git clone --recursive https://github.com/jts/nanopolish.git \
    && pip3 install -r nanopolish/scripts/requirements.txt \
    && cd nanopolish && make --silent --ignore-errors && cd .. \
# Clair & Clair-env
    && wget -q https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && bash Miniconda3-latest-Linux-x86_64.sh -b -p /opt/conda \
    && rm Miniconda3-latest-Linux-x86_64.sh
RUN conda config --add channels bioconda \
    && conda config --add channels conda-forge \
    && conda create -n clair-env -c bioconda -y \
    clair tabix r-sys bioconductor-dss \
    && conda clean --all
RUN echo "source activate clair-env" >> /etc/profile
ENV PATH=/opt/conda/envs/clair-env/bin:$PATH
RUN /bin/bash -c ". activate clair-env \
    && pypy3 -m ensurepip \
    && pip3 install --upgrade pip setuptools wheel \
    && pypy3 -m pip install --no-cache-dir intervaltree \
    && conda deactivate \
    && mkdir ont && cd ont \
    && wget -q http://www.bio8.cs.hku.hk/clair_models/ont/122HD34.tar \
    && tar -xf 122HD34.tar && rm 122HD34.tar && cd .. \
    && chmod -R +x /opt \
    && ln -rs /opt/conda/envs/clair-env/bin/clair.py /bin/clair \
    && ln -rs /opt/conda/envs/clair-env/bin/tabix /bin/tabix \
    && ln -rs /opt/clair/nanopolish/nanopolish /bin/nanopolish"
ENTRYPOINT ["/opt/clair/entrypoint.sh"]
