FROM nfcore/base

MAINTAINER Diego Brambilla
LABEL description="Docker image containing all requirements for nf-core/ampliseq pipeline, including dependecies for R-DADA2 scripts"

## Install necessary tools
COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a && \
    cd /opt && git clone https://github.com/erikrikarddaniel/eemisdada2.git && \
    ln -s /opt/eemisdada2/src/R/dada2* /usr/local/bin
ENV PATH /opt/conda/envs/nf-core-ampliseq-1.1.0/bin:$PATH
## Required to build the container properly
RUN mkdir -p /root/.config/matplotlib
RUN echo "backend : Agg" > /root/.config/matplotlib/matplotlibrc
## Don't recache on each execution, do that once per build process
RUN qiime dev refresh-cache

