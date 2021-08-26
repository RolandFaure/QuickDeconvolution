# QuickDeconvolution

Quick and scalable software to deconvolve read clouds from linked-reads experiments. When several fragments of DNA have been sequenced with the same barcode, QuickDeconvolution provides the user with enhanced barcodes to distinguish the reads coming from the different fragments

## Installation

You will need make and cmake to compile the sources. In the desired folder, run

```
git clone https://github.com/RolandFaure/QuickDeconvolution.git
cd QuickDeconvolution/
cmake ./
make
```

An executable named QuickDeconvolution should appear in the folder.

## Usage

```bash
SYNOPSIS
        ./QuickDeconvolution -i [<input-file>] -o [<output-file>] [-k [<k>]] [-w
            [<w>]] [-d [<d>]] [-t [<t>]] [-a [<a>]]

OPTIONS
        -k, --kmers-length
                    size of kmers [default:20]

        -w, --window-size
                    size of window guaranteed to contain at least one minimizing kmer [default:40]

        -d, --density
                    on average 1/2^d kmers are sparse kmers [default:3]

        -t, --threads
                    number of threads [default:1]

        -a, --dropout
                    QD does not try to deconvolve clouds smaller than this value [default:0]

```

### Input

QuickDeconvolution takes as input `-i` a fasta or a fastq file containing barcoded reads with the tag `BX:Z` designating a barcode (this is the default output of [longranger basic](https://support.10xgenomics.com/genome-exome/software/pipelines/latest/advanced/other-pipelines)). For example
```
@read_456 cov:23.45 BX:Z:AAAACTGTAT
```
If the reads are paired, provide QuickDeconvolution with an interleaved file where the two ends of the pairs have the same name, it will recognize it.

### Output

QuickDeconvolution outputs the fasta/q file given as input, with an additional tag (-0, -1, -2...) at the end of the line, so that the deconvolved reads look like 
```
@read_456 cov:23.45 BX:Z:AAAACTGTAT-1
```
Within each barcode, reads having the same tag come from the same fragment. WARNING: the -0 tag is a special tag, indicating reads that could not be deconvolved by the program.
