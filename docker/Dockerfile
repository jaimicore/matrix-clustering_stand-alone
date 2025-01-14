FROM r-base:4.4.2

## --------------------------- ##
## Install system requirements ##
## --------------------------- ##

RUN mkdir -p /opt/software
WORKDIR /opt

RUN apt-get update -qq && \
        apt-get install -y \
                apt-utils=2.9.21

RUN apt-cache policy \
                apt-utils \
                build-essential \
                freetype2-demos \
                freetype2-doc \
                git \
                libcurl4-openssl-dev \
                libcairo2-dev \
                libfontconfig1-dev \
                libfreetype-dev \
                libfreetype6 \
                libfreetype6-dev \
                libmysqlclient-dev \
                libnss-unknown \
                libpq-dev \
                librsvg2-bin \
                librsvg2-dev \
                libudunits2-dev \
                littler \
                make \
                pkg-config \
                procps \
                libsharpyuv-dev \
                r-cran-svglite \
        r-cran-systemfonts

                
RUN apt-get update -qq && \
        apt-get install -y --allow-downgrades \
                apt-utils=2.9.21 \
                build-essential=12.12 \
                libwebp7=1.4.0-0.1+b1 \
                libsharpyuv0=1.4.0-0.1+b1 \
                libsharpyuv-dev=1.4.0-0.1+b1 \
                freetype2-demos=2.13.3+dfsg-1 \
                freetype2-doc=2.13.3+dfsg-1 \
                git=1:2.45.2-1 \
                libcurl4-openssl-dev=8.11.1-1 \
                libcairo2-dev=1.18.2-2 \
                libfontconfig1-dev=2.15.0-2 \
                libfreetype-dev=2.13.3+dfsg-1 \
                libfreetype6=2.13.3+dfsg-1 \
                # libfreetype6-dev \
                libmysqlclient-dev=8.0.40-2 \
                libnss-unknown=0.0.2-4+b2 \
                libpq-dev=17.2-1 \
                librsvg2-bin=2.59.1+dfsg-1 \
                librsvg2-dev=2.59.1+dfsg-1 \
                libudunits2-dev=2.2.28-7 \
                littler=0.3.20-2 \
                make=4.4.1-1 \
                pkg-config=1.8.1-4 \
                procps=2:4.0.4-6


RUN apt-get install --yes --quiet \
        r-cran-svglite=2.1.3-1+b1 \
        r-cran-systemfonts=1.1.0-1

## ------------------- ##
## Copy required files ##
## ------------------- ##

RUN git clone https://github.com/jaimicore/matrix-clustering_stand-alone.git
COPY docker/requirements.R /opt/requirements.R

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