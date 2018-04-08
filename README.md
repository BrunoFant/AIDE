AID: annotation-based isoform discovery and quantification from RNA-seq data
================
Wei Vivian Li, Jingyi Jessica Li
2018-04-08

<!-- README.md is generated from README.Rmd. Please edit that file -->
Latest News
-----------

> 2018/04/08:

-   Version 0.0.1 is released!

Introduction
------------

AID is a statistical method which identifies full-length mRNA isoforms from a novel perspective: using the likelihood ratio test to find novel isoforms in a stepwise manner given annotated isoforms, by prioritizing and selectively borrowing information from the annotated isoforms.

Any suggestions on the package are welcome! For technical problems, please report to [Issues](https://github.com/Vivianstats/AID/issues). For suggestions and comments on the method, please contact Wei (<liw@ucla.edu>) or Dr. Jessica Li (<jli@stat.ucla.edu>).

Installation
------------

You can install `AID` from github with:

``` r
# install.packages("devtools")
devtools::install_github("Vivianstats/AID")
```

Quick start
-----------

`aid` requires three input files:

-   The GTF file of the genome annotation;
-   The BAM file of the RNA-seq sample. The BAM file should be sorted and the index BAI file should be present in the same directory as the BAM file;
-   The FASTA file of the genome sequences.

The final output of `aid` is a GTF file named "transcripts.gtf", which contains the reconstructed transcripts and their corresponding abudance levels. The package has been tested using the [GENCODE annotation](https://www.gencodegenes.org/releases/24.html). This is a basic example which shows how to use the `aid` function.

``` r
aid(gtf_path = "./hg19.gtf",     #full path of the GTF file
    bam_path = "./example.bam",  #full path of the BAM file
    fasta_path = "./hg19.fa",    #full path of the FASTA file
    out_dir = "./",              #output directory of temporary and filnal results
    readLen = 100,               #read length used to calculate ioform effective length
    strandmode = 0,              #library type of the RNA-seq sample
    ncores = 20                  #number of cores used for parallel computation 
    )
```

Please refer to the package [manual](https://github.com/Vivianstats/AID/blob/master/inst/docs/) for a full list of arguments and detailed usage.
