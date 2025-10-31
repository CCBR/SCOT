FROM nciccbr/ccbr_ubuntu_22.04:v4

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

# Pin channels and update
RUN conda config --add channels conda-forge \
 && conda config --add channels bioconda \
 && conda config --set channel_priority strict

# install conda packages
RUN mamba install -y -c conda-forge \
    r-base=${R_VERSION} \
    r-devtools \
  && conda clean -afy

# install R package
COPY . /opt2/SCOT
RUN R -e "devtools::install_local('/opt2/SCOT', dependencies = TRUE, repos='http://cran.rstudio.com')"

# Save Dockerfile in the docker
COPY Dockerfile /opt2/Dockerfile_${REPONAME}.${BUILD_TAG}
RUN chmod a+r /opt2/Dockerfile_${REPONAME}.${BUILD_TAG}

# cleanup
WORKDIR /data2
RUN apt-get clean && apt-get purge \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
