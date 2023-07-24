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

WORKDIR /opt
RUN git clone https://github.com/jaimicore/matrix-clustering_stand-alone.git
COPY requirements.R /opt/requirements.R

## ---------------------------- ##
## Install package requirements ##
## ---------------------------- ##

RUN Rscript requirements.R
ENV R_LIBS=${R_LIBS}:/opt/software

RUN cd matrix-clustering_stand-alone/compare-matrices-quick \
	&& make \
	&& cd ../..

## --------------------------------------------------- ##
## Create symlinks to required scripts and update PATH ##
## --------------------------------------------------- ##

RUN chmod +x /opt/matrix-clustering_stand-alone/matrix-clustering.R \
    && ln -s /opt/matrix-clustering_stand-alone/matrix-clustering.R /usr/bin/matrix-clustering

## -------------------------- ##
## Switch to non-root account ##
## -------------------------- ##

RUN useradd -u 4444 user \
    && chown -R user /opt/ \
    && chgrp -R user /opt/

USER user