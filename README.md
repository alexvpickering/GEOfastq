# GEOfastq

### Install GEOfastq

To download and install `GEOfastq`:

```R
install.packages('remotes')
remotes::install_github('alexvpickering/GEOfastq')
```

### Install Aspera Connect

`GEOfastq` uses [aspera connect](https://downloads.asperasoft.com/en/downloads/8?list) because it is faster than ftp. Download and install it according to the [documentation](https://downloads.asperasoft.com/en/documentation/8). For me (Fedora 30), this works:

```bash
wget https://download.asperasoft.com/download/sw/connect/3.9.6/ibm-aspera-connect-3.9.6.173386-linux-g2.12-64.tar.gz
tar -zxvf ibm-aspera-connect-3.9.6.173386-linux-g2.12-64.tar.gz
./ibm-aspera-connect-3.9.6.173386-linux-g2.12-64.sh
```

I also had to make sure `ascp` was on the the `PATH`:

```bash
echo 'export PATH=$HOME/.aspera/connect/bin:$PATH' >> ~/.bashrc
source ~/.bashrc
```

For Rstudio to find `ascp` on the `PATH`, I also had to add this to a .Renviron:

```bash
echo 'PATH=${HOME}/.aspera/connect/bin:${PATH}' >> ./Renviron
```

After restarting Rstudio, to confirm things are set up properly:


```R
# should have the above path added
Sys.getenv('PATH')

# should print info about Aspera Connect
system('ascp --version')
```

### Install docker image

To install `GEOfastq` and Aspera Connect from a pre-built docker image:

```bash
# retrieve pre-built geofastq docker image
wget https://drugseqr.s3.us-east-2.amazonaws.com/geofastq_latest.tar.gz
sudo docker load < geofastq_latest.tar.gz

# run interactive container with host portion of `-v host:container` mounted where you want to persist data to
sudo docker run -it --rm \
  -v /srv/shiny-server:/srv/shiny-server \
  geofastq /bin/bash
```


### Usage

First crawl a study page on [GEO](https://www.ncbi.nlm.nih.gov/geo/) to get study metadata and corresponding fastq.gz download links on [ENA](https://www.ebi.ac.uk/ena):

```R
gse_name <- 'GSE117570'
srp_meta <- GEOfastq::get_srp_meta(gse_name)
```

Next, subset `srp_meta` to samples that you want, then download:

```R
srp_meta <- srp_meta[srp_meta$source_name == 'Adjacent normal', ]
GEOfastq::get_fastqs(gse_name, srp_meta)
```

That's all folks!

