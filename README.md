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

An executable named QuickDeconvolution should appear in the folder. A small test file "test.fastq", from a simulated sequencing experiment on a small synthetic genome, is provided to run the program.

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
                    on average 1/2^d kmers are indexed [default:3]

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

### Options

Option -a is the dropout option: the program disregard all clouds containing fewer reads than this value. You may want to use the option if you know you'll need clouds of a certain size for your downstream analyses, in which case it might be a waste of time to deconvolve the smallest clouds.

Option -t is the number of threads to launch simultaneously on the program. Wall-clock time decreases and RAM usage increases with the number of threads.

Options k, w and d are parameters of the alignment within QuickDeconvolution. The deconvolution should not be very sensitive to these values. 
k is the length of the k-mers. Avoid decreasing k below 15. 
d is to monitor the density of sparse k-mers. On average 1/2^d k-mers will be sparse.
While choosing sparse k-mers, the program is ensured to choose at least 1 k-mer in a window of size w.
Decreasing w and/or d will increase precision at the expense of run-time. Keep w within the range [10,50] and d within range [1,5].
