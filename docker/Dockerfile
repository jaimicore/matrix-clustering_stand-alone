FROM r-base:4.4.2

## --------------------------- ##
## Install system requirements ##
## --------------------------- ##

RUN mkdir -p /opt/software
WORKDIR /opt

RUN apt-get update -qq && \
        apt-get install -y \
                apt-utils \
                build-essential \
                git \
                littler \
                make \
                libcurl4-openssl-dev \
                pkg-config \
                libfontconfig1-dev \
                libfreetype-dev \
                libfreetype6 \
                libfreetype6-dev \
                libcairo2-dev \
                freetype* \
                libudunits2-dev \
                librsvg2-dev \
                librsvg2-bin \
                libpq-dev \
                libmysqlclient-dev \
                libnss-unknown \
                procps

RUN apt-get install --yes --quiet \
        r-cran-svglite \
        r-cran-systemfonts

## ------------------- ##
## Copy required files ##
## ------------------- ##

RUN git clone https://github.com/jaimicore/matrix-clustering_stand-alone.git
COPY requirements.R /opt/requirements.R

## ---------------------------- ##
## Install package requirements ##
## ---------------------------- ##

RUN Rscript requirements.R

WORKDIR /opt/matrix-clustering_stand-alone/compare-matrices-quick
RUN make
WORKDIR /opt/software/matrix-clustering_stand-alone

RUN fc-cache -rvf && rm -rf /var/cache/*
ENV OSFONTDIR=/usr/share/texmf/fonts/

## --------------------------------------------------- ##
## Create symlinks to required scripts and update PATH ##
## --------------------------------------------------- ##

RUN chmod +x /opt/matrix-clustering_stand-alone/matrix-clustering.R \
        && ln -s /opt/matrix-clustering_stand-alone/matrix-clustering.R /usr/bin/matrix-clustering

RUN chmod +x /opt/matrix-clustering_stand-alone/convert-matrix.R \
        && ln -s /opt/matrix-clustering_stand-alone/convert-matrix.R /usr/bin/convert-matrix

## -------------------------- ##
## Switch to non-root account ##
## -------------------------- ##

ENV HOME=/home

CMD ["/bin/bash"]