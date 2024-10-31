FROM fedora:37

## --------------------------- ##
## Install system requirements ##
## --------------------------- ##

RUN mkdir -p /opt/software
COPY requirements-sys.txt /opt/requirements-sys.txt
WORKDIR /opt
RUN xargs dnf -y install < "requirements-sys.txt"

## ------------------- ##
## Copy required files ##
## ------------------- ##

RUN git clone https://github.com/jaimicore/matrix-clustering_stand-alone.git
COPY requirements.R /opt/requirements.R

## ---------------------------- ##
## Install package requirements ##
## ---------------------------- ##

COPY data.table_1.14.8.tar.gz /opt/data.table_1.14.8.tar.gz

RUN Rscript requirements.R
ENV R_LIBS=${R_LIBS}:/opt/software

WORKDIR /opt/matrix-clustering_stand-alone/compare-matrices-quick
RUN make
WORKDIR /opt/software/matrix-clustering_stand-alone

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

RUN useradd -u 4444 user \
    && chown -R user /opt/ \
    && chgrp -R user /opt/

USER user