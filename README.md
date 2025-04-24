
# RVcov

<!-- badges: start -->
<!-- badges: end -->

The goal of RVcov is to perform rare-variant aggregation testing when users have
RAREMETALWORKER formatted summary statistics and a VCF containing a set of individuals
to serve as an external reference panel. The agg_test() function will estimate 
the covariance of the score statistics using the VCF from the external reference 
panel and can perform various rare-variant aggregation tests using this estimated 
covariance. Options to supply either VEP annotation files or user-made RAREMETAL 
formatted group files. Option also to perform the two-stage approach whereby covariance
is estimated only when using the null covariance gives a -log10(p-value) above some 
threshold (default to 3). 

A tutorial of using the package can be found in the vignettes. After package installation,
a copy of the data 
used in the vignette as well as another README containing a description of how the 
data was generated is available in the extdata/ folder. The names of the available 
files can be found in the Tutorial vignette, and a description on how to
access these files is available at the following link: 

https://r-pkgs.org/data.html#sec-data-extdata

Please note that this package works only on Linux, and you must have previously 
installed tabix, RAREMETAL, and BCFtools. The tabix and BCFtools must be able to 
be called using the commands tabix and bcftools respectively. By default we assume
RAREMETAL can be called using the command raremetal but there is an option to supply
an alternate path/name to call raremetal.

A detailed explanation of our method is available in our paper: (insert paper link eventually)

## Installation

Note that our package requires the VariantAnnotation package to be installed 
and is not installed automatically. You can install both VariantAnnotation and 
the development version of RVcov from [GitHub](https://github.com/balr411) with:

``` r
# install.packages("devtools")
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("VariantAnnotation")
devtools::install_github("balr411/RVcov")
```

Instructions on how to install tabix (from SAMtools) and BCFtools can be found at the following link:

http://www.htslib.org/download/

We recommend installing the following alternative version of RAREMETAL as there are
multiple bugs with the traditional version:

https://github.com/statgen/nullmetal

Instructions on installation are also available at that link.




