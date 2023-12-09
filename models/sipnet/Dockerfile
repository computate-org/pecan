# this needs to be at the top, what version are we building
ARG IMAGE_VERSION="latest"
ARG FROM_ORG="pecan"
ARG FROM_IMAGE="base"

# ----------------------------------------------------------------------
# BUILD SIPNET BINARY
# ----------------------------------------------------------------------
FROM ${FROM_ORG}/${FROM_IMAGE}:${IMAGE_VERSION} as model-binary

# Some variables that can be used to set control the docker build
ARG MODEL_VERSION=git
ARG GITHUB_ORG="PecanProject"
ARG GITHUB_REPO="sipnet"

# download, unzip and build sipnet
WORKDIR /src
RUN echo "https://github.com/${GITHUB_ORG}/${GITHUB_REPO}.git"
RUN git clone "https://github.com/${GITHUB_ORG}/${GITHUB_REPO}.git" \
    && cd sipnet \
    && if [ "${MODEL_VERSION}" != "git" ]; then git checkout ${MODEL_VERSION}; fi \
    && env LDLAGS=-lm make

########################################################################

# ----------------------------------------------------------------------
# BUILD PECAN FOR SIPNET
# ----------------------------------------------------------------------
FROM ${FROM_ORG}/${FROM_IMAGE}:${IMAGE_VERSION}

# ----------------------------------------------------------------------
# INSTALL SIPNET SPECIFIC PIECES
# ----------------------------------------------------------------------

# INSTALL PEcAn.data.atmosphere
#COPY modules/data.atmosphere /src/data.atmosphere/
#RUN Rscript --vanilla -e  "devtools::install_local('/src/data.atmosphere'); q(status=as.integer(!'PEcAn.data.atmosphere' %in% installed.packages()))" 

# INSTALL PEcAn.SIPNET
#COPY models/sipnet /src/sipnet/
#RUN Rscript --vanilla -e  "devtools::install_local('/src/sipnet'); q(status=as.integer(!'PEcAn.SIPNET' %in% installed.packages()))" 

# ----------------------------------------------------------------------
# SETUP FOR SPECIFIC SIPNET VERSION
# ----------------------------------------------------------------------

# Some variables that can be used to set control the docker build
ARG MODEL_VERSION=git

# Setup model_info file
COPY model_info.json /work/model.json
RUN sed -i -e "s/@VERSION@/${MODEL_VERSION}/g" \
           -e "s#@BINARY@#/usr/local/bin/sipnet.${MODEL_VERSION}#g" /work/model.json

# COPY sipnet binary
COPY --from=model-binary /src/sipnet/sipnet /usr/local/bin/sipnet.${MODEL_VERSION}
