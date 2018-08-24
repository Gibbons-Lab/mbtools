FROM rocker/rstudio
LABEL author="Christian Diener <mail[at]cdiener.com>"

RUN sudo apt-get update && sudo apt-get install -y --no-install-recommends \
    zlib1g-dev bowtie2 samtools

# Setup dependencies
RUN Rscript -e "source('http://bioconductor.org/biocLite.R'); \
    biocLite('BiocInstaller'); setRepositories(ind=1:2); \
    install.packages('devtools'); \
    devtools::install_github('Gibbons-Lab/mbtools')" \
    && rm -rf /tmp/*

RUN mkdir /data
COPY docs/mock_example.Rmd /data
RUN chown -R rstudio:rstudio /data
