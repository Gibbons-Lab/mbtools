FROM rocker/rstudio
LABEL author="Christian Diener <mail[at]cdiener.com>"

RUN sudo apt-get update && sudo apt-get install -y --no-install-recommends \
    zlib1g-dev bowtie2 samtools $$ \
    curl -L https://github.com/lh3/minimap2/releases/download/v2.13/minimap2-2.13_x64-linux.tar.bz2 | tar -jxvf - /usr/bin/minimap2

# Setup dependencies
RUN Rscript -e "source('http://bioconductor.org/biocLite.R'); \
    biocLite('BiocInstaller'); setRepositories(ind=1:2); \
    install.packages('devtools'); \
    devtools::install_github('Gibbons-Lab/mbtools')" \
    && rm -rf /tmp/*

RUN mkdir /data
COPY docs/mock_example.Rmd /data
RUN chown -R rstudio:rstudio /data
