FROM rocker/r-ubuntu:20.04
MAINTAINER vladimir.kovacevic@wonderland.ai
RUN apt-get update &&  apt install -y  git-all

RUN R -e "install.packages('argparser', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('huge', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('MASS', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('rjson', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('stringr', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('fs', repos='http://cran.rstudio.com/')"
#RUN R -e "install.packages('CTD', repos='http://cran.rstudio.com/')"

WORKDIR /opt
# Avoid pulling from docker cache in order to get the latest code from GitHub
ADD "https://www.random.org/cgi-bin/randbyte?nbytes=10&format=h" skipcache
RUN git clone https://github.com/BRL-BCM/CTD.git 
WORKDIR /opt/CTD

COPY Dockerfile /opt
