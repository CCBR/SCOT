FROM nciccbr/ccbr_ubuntu_base_20.04:v7

# build time variables
ARG BUILD_DATE="000000"
ENV BUILD_DATE=${BUILD_DATE}
ARG BUILD_TAG="000000"
ENV BUILD_TAG=${BUILD_TAG}
ARG REPONAME="000000"
ENV REPONAME=${REPONAME}

ARG R_VERSION=4.3.2
ENV R_VERSION=${R_VERSION}

SHELL ["/bin/bash", "-lc"]

# install system dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    libcairo2-dev \
    libcurl4-openssl-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libfribidi-dev \
    libharfbuzz-dev \
    libhdf5-dev \
    libssl-dev \
    libsz2 \
    libxml2-dev \
    pkg-config

# install conda packages
RUN mamba install -y -c conda-forge -c bioconda \
    r-base=${R_VERSION} \
    r-devtools \
    bioconductor-annotationdbi \
    bioconductor-aucell \
    bioconductor-celldex \
    bioconductor-dropletutils \
    bioconductor-gseabase \
    bioconductor-ontoProc \
    bioconductor-org.hs.eg.db \
    bioconductor-org.mm.eg.db \
    bioconductor-orthology.eg.db \
    bioconductor-scdblfinder \
    bioconductor-singlecellexperiment \
    bioconductor-singleR \
    r-cffr \
    r-clustersim \
    r-dplyr \
    r-ggplot2 \
    r-ggpubr \
    r-ggrepel \
    r-goodpractice \
    r-harmony \
    r-knitr \
    r-lintr \
    r-magrittr \
    r-rcmdcheck \
    r-reshape2 \
    r-rlang \
    r-rmarkdown \
    r-roxygen2 \
    r-r.utils \
    r-seurat \
    r-seuratobject \
    r-testthat \
    r-usethis \
  && conda clean -afy

# install R package
COPY . /opt2/SCOT
RUN R -e "devtools::install_local('/opt2/SCOT', dependencies = TRUE, repos='http://cran.rstudio.com')" && \
  R -e "library(SCOT)" && \
  R -s -e "readr::write_csv(tibble::as_tibble(installed.packages()), '/mnt/r-packages.csv')"

# Save Dockerfile in the docker
COPY Dockerfile /opt2/Dockerfile_${REPONAME}.${BUILD_TAG}
RUN chmod a+r /opt2/Dockerfile_${REPONAME}.${BUILD_TAG}

# cleanup
WORKDIR /data2
RUN apt-get clean && apt-get purge \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
