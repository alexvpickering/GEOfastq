# GEOfastq

### Install GEOfastq

To download and install `GEOfastq`:

```R
install.packages('remotes')
remotes::install_github('alexvpickering/GEOfastq')
```

### Install Aspera Connect (optional)

`GEOfastq` can use [aspera
connect](https://downloads.asperasoft.com/en/downloads/8?list) to download
fastqs. It is faster than ftp for large single-file downloads (single-cell
fastqs).
To download and install it according to the
[documentation](https://downloads.asperasoft.com/en/documentation/8). For me
(Fedora 30), this works:

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
system2('ascp', '--version')
```

### Install docker image

To install `GEOfastq` and Aspera Connect from a pre-built docker image:

```bash
# retrieve pre-built geofastq docker image
docker pull alexvpickering/geofastq

# run interactive container with host portion of 
#`-v host:container` mounted where you want to persist data to
sudo docker run -it --rm \
  -v /srv:/srv \
  geofastq /bin/bash
```


### Usage

First crawl a study page on [GEO](https://www.ncbi.nlm.nih.gov/geo/) to get
study metadata and corresponding fastq.gz download links on
[ENA](https://www.ebi.ac.uk/ena):

```R
library(GEOfastq)

gse_name <- 'GSE117570'
#' gse_text <- crawl_gse(gse_name)
#' gsm_names <- extract_gsms(gse_text)
#' srp_meta <- crawl_gsms(gsm_names)
```

Next, subset `srp_meta` to samples that you want, then download:

```R
srp_meta <- srp_meta[srp_meta$source_name == 'Adjacent normal', ]
get_fastqs(srp_meta, data_dir = tempdir())
```

That's all folks! GOTO: `kallisto`?

