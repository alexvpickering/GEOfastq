# GEOfastq

## Installation

To download and install `GEOfastq`:

```R
remotes::install_github('alexvpickering/GEOfastq')
```

`GEOfastq` uses [aspera connect](https://downloads.asperasoft.com/en/downloads/8?list). Download and install it according to the [documentation](https://downloads.asperasoft.com/en/documentation/8). For me (Fedora 30), this works:

```bash
wget https://download.asperasoft.com/download/sw/connect/3.9.6/ibm-aspera-connect-3.9.6.173386-linux-g2.12-64.tar.gz
tar -zxvf ibm-aspera-connect-version.tar.gz ibm-aspera-connect-3.9.6.173386-linux-g2.12-64.tar.gz
./ibm-aspera-connect-3.9.6.173386-linux-g2.12-64.sh
```

I also had to make sure `ascp` was on the the `PATH`:

```bash
echo 'export PATH=$HOME/.aspera/connect/bin/:$PATH' >> ~/.bashrc
source ~/.bashrc
```

For Rstudio to find `ascp` on the `PATH`, I also had to add this to a .Renviron:

```bash
echo 'PATH=${HOME}/.aspera/connect/bin/:${PATH}' >> ./Renviron
```

After restarting Rstudio, to confirm things are set up properly:


```R
# should have the above path added
Sys.getenv('PATH')

# should print info about Aspera Connect
system('ascp --version')
```


## Installation

## Usage

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

