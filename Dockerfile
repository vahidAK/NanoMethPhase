FROM ubuntu:focal
# NanoMethPhase & SNVoter
ARG DEBIAN_FRONTEND=noninteractive
RUN apt update && apt install -y wget gcc git libz-dev build-essential \
    dirmngr apt-transport-https ca-certificates software-properties-common \
    python3 python3-pip make python-dev libhdf5-dev libcudart10.1 \
    && pip3 install --upgrade pip setuptools wheel \
    && pip3 install nanomethphase \
    && apt install -y --no-install-recommends r-base r-cran-biocmanager \
    && apt-get clean \
    && Rscript -e "BiocManager::install('DSS')" \
# WhatsHap
    && pip3 install whatshap \
# non-root user
    && groupadd --gid 5000 newuser \
    && useradd --home-dir /home/newuser --create-home --uid 5000 \
        --gid 5000 --shell /bin/sh --skel /dev/null newuser
USER newuser
# NanoPolish
RUN cd /home/newuser \
    && git clone --recursive https://github.com/jts/nanopolish.git \
    && pip3 install -r nanopolish/scripts/requirements.txt \
    && cd nanopolish && make && cd .. \
# Clair
    && cd /home/newuser \
    && wget -q https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && bash Miniconda3-latest-Linux-x86_64.sh -b \
    && rm Miniconda3-latest-Linux-x86_64.sh
ENV PATH="/home/newuser/miniconda3/bin:${PATH}"
ARG PATH="/home/newuser/miniconda3/bin:${PATH}"
RUN conda install -c conda-forge -c bioconda -y python=3.7 clair tabix \
    && /home/newuser/miniconda3/bin/pypy3 -m ensurepip \
    && /home/newuser/miniconda3/bin/pypy3 -m pip install \
        --no-cache-dir intervaltree \
    && conda init bash \
    && echo ". /home/newuser/miniconda3/etc/profile.d/conda.sh" >> /home/newuser/.bashrc \
    && echo "conda activate base" >> /home/newuser/.bashrc \
    && mkdir /home/newuser/ont && cd /home/newuser/ont \
    && wget -q http://www.bio8.cs.hku.hk/clair_models/ont/122HD34.tar \
    && tar -xf 122HD34.tar && rm 122HD34.tar && cd ../
# NanoPolish & Clair in path
USER root
RUN chmod a+x /home/newuser/miniconda3/bin/clair.py \
    && ln -rs /home/newuser/miniconda3/bin/clair.py /bin/clair \
    && chmod a+x /home/newuser/miniconda3/bin/tabix \
    && ln -rs /home/newuser/miniconda3/bin/tabix /bin/tabix \
    && chmod a+x /home/newuser/nanopolish/nanopolish \
    && ln -rs /home/newuser/nanopolish/nanopolish /bin/nanopolish
ENTRYPOINT . /home/newuser/miniconda3/etc/profile.d/conda.sh
