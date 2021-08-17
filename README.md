# QuickDeconvolution

Quick and scalable software to deconvolve read clouds from linked-reads experiments. When sevelal fragments of DNA have been sequenced with the same barcode, QuickDeconvolution provides the user with enhanced barcodes to distinguish the reads coming from the different fragments

## Usage

```bash
SYNOPSIS
        ./QuickDeconvolution -i [<i>] -f [<f>] -o [<o>] [-k [<k>]] [-w
            [<w>]] [-d [<d>]] [-t [<t>]] [-a [<a>]]

OPTIONS
        -k, --kmers-length
                    size of kmers

        -w, --window-size
                    size of window guaranteed to contain at least one minimizing kmer

        -d, --density
                    on average 1/2^d kmers are sparse kmers

        -t, --threads
                    number of threads

        -a, --dropout
                    QD does not try to deconvolve clouds smaller than this value [default:0]

```

### Input

QuickDeconvolution takes as input a fasta or a fastq file containing barcoded reads with the tag `BX:Z` designating a barcode (this is the default output of [longranger basic](https://support.10xgenomics.com/genome-exome/software/pipelines/latest/advanced/other-pipelines)). For example
```
@read_456 cov:23.45 BX:Z:AAAACTGTAT
```
If the reads are paired, provide QuickDeconvolution with an interleaved file where the two ends of the pairs have the same name, it will recognize it.
