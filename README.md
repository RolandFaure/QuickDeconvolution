# QuickDeconvolution

Quick and scalable software to deconvolve read clouds from linked-reads experiments without a reference genome. When several fragments of DNA have been sequenced with the same barcode, QuickDeconvolution provides the user with enhanced barcodes to distinguish the reads coming from the different fragments

## Installation

You can install QuickDeconvolution through Bioconda
```
conda install -c bioconda quickdeconvolution
```

## Build from source

QuickDeconvolution requires CMake 3.20 or newer, a C++11-capable compiler,
and Boost.Iostreams with gzip support. Install the build dependencies for
your platform before running CMake.
<details>
<summary>
    Installing boost
</summary>

---

**Linux**
```bash
# Ubuntu or Debian-derived systems
sudo apt update
sudo apt install build-essential cmake libboost-iostreams-dev

# Arch Linux
sudo pacman -S base-devel cmake boost
```

**macOS with Homebrew**

Install the Xcode Command Line Tools with `xcode-select --install` if you
haven’t already, then run:

```bash
brew install cmake boost
```
---
</details>

### Compile

```bash
git clone https://github.com/BLinkseq/QuickDeconvolution.git
cd QuickDeconvolution
cmake -S . -B build
cmake --build build --parallel
```

An executable named QuickDeconvolution should appear in the folder. A small test file "test.fastq.gz", from a
simulated sequencing experiment on a small synthetic genome, is provided in the folder "test_data" to test the program.
```
QuickDeconvolution -i test_data/test.fastq -o test_data/test_out.fastq
```
The program should run in less than a minute and output in test_out.fastq the reads, with barcode extensions (-1, -2,...). This is only intended as a test to see if QD is running: the deconvolution is expected to be very bad because the synthetic genome is very short (thus two long reads overlap with high probability).

## Usage

```bash
SYNOPSIS
    ./QuickDeconvolution -i [<input-file>] -o [<output-file>] [-k [<k>]] [-w
        [<w>]] [-d [<d>]] [-t [<t>]] [-a [<a>]]

OPTIONS
    -k, --kmers-length
        size of kmers [default:21]
    
    -w, --window-size
        size of window guaranteed to contain at least one minimizing kmer [default:40]
    
    -d, --density
        on average 1/2^d kmers are indexed [default:3]
    
    -t, --threads
        number of threads [default:1]
    
    -a, --dropout
        QD does not try to deconvolve clouds smaller than this value [default:0]
    -m, --metagenome
        Use this option on metagenomic samples
    
    -v, --verbose
        Log run information to stderr
```

### Input

QuickDeconvolution takes as input `-i` a fasta or a fastq file containing barcoded reads with the tag `BX:Z` designating a barcode (this is the default output of [longranger basic](https://support.10xgenomics.com/genome-exome/software/pipelines/latest/advanced/other-pipelines)). For example
```
@read_456 cov:23.45 BX:Z:AAAACTGTAT
```
If the reads are paired, provide QuickDeconvolution with an interleaved file where the two ends of the pairs have the same name, it will recognize it. To interleave two files, you can use this command line: 
```
paste -d '\n' <(awk '{if (NR%4==1){printf"\n";printf $0;} else{printf "((()))"$0;}}' reads_foward.fq) <(awk '{if (NR%4==1){printf"\n";printf $0;}else{ printf "((()))"$0;}}' reads_reverse.fq) | sed 's/((()))/\n/g' > sequencing_reads_interleaved.fastq
```

### Output

QuickDeconvolution outputs the fasta/q file given as input, with an additional tag (-0, -1, -2...) at the end of the line, so that the deconvolved reads look like 
```
@read_456 cov:23.45 BX:Z:AAAACTGTAT-1
```
Within each barcode, reads having the same tag come from the same fragment. WARNING: the -0 tag is a special tag, indicating reads that could not be deconvolved by the program. If a tag is already present, QuickDeconvolution will nonetheless append a new tag at the end of the barcode:
```
@read_456 cov:23.45 BX:Z:AAAACTGTAT-1-3
```

### Options

Option -a is the dropout option: the program disregard all clouds containing fewer reads than this value. You may want to use the option if you know you'll need clouds of a certain size for your downstream analyses, in which case it might be a waste of time to deconvolve the smallest clouds.

Option -t is the number of threads to launch simultaneously on the program. Wall-clock time decreases and RAM usage increases with the number of threads.

Options k, w and d are parameters of the alignment within QuickDeconvolution. The deconvolution should not be very sensitive to these values. 
k is the length of the k-mers. Avoid decreasing k below 15. 
d is to monitor the density of sparse k-mers. On average 1/2^d k-mers will be sparse.
While choosing sparse k-mers, the program is ensured to choose at least 1 k-mer in a window of size w.
Decreasing w and/or d may in some cases increase precision at the expense of run-time. Keep w within the range [10,50] and d within range [1,5].

## License

QuickDeconvolution is distributed under the license GPL3

## Citation

QuickDeconvolution is published in *Bioinformatics advances*. You can cite using:
Faure, Roland, and Dominique Lavenier. “QuickDeconvolution: Fast and Scalable Deconvolution of Linked-Read Sequencing Data.” *Bioinformatics Advances*, September 26, 2022, vbac068. _https://doi.org/10.1093/bioadv/vbac068_.
