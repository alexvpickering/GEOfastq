FROM rocker/r-ver:3.6.1

RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    wget && \
    rm -rf /var/lib/apt/lists/*


# install aspera connect and GEOfastq
RUN wget https://download.asperasoft.com/download/sw/connect/3.9.6/ibm-aspera-connect-3.9.6.173386-linux-g2.12-64.tar.gz && \
    tar -zxvf ibm-aspera-connect-3.9.6.173386-linux-g2.12-64.tar.gz && \
    ./ibm-aspera-connect-3.9.6.173386-linux-g2.12-64.sh

RUN R -e "install.packages('remotes')" && \
    R -e "remotes::install_github('rstudio/renv@0.7.1-20')"

COPY renv.lock .

RUN R -e 'options(renv.consent = TRUE); renv::restore()' && \
    R -e "remotes::install_github('alexvpickering/GEOfastq', dependencies = FALSE, upgrade = FALSE)"

ENV PATH="/root/.aspera/connect/bin:$PATH"

# save image to a tar.gz file and upload to s3
# sudo docker build -t geofastq:latest .
# sudo docker save geofastq:latest | gzip > geofastq_latest.tar.gz
# aws s3 cp geofastq_latest.tar.gz s3://drugseqr/geofastq_latest.tar.gz
