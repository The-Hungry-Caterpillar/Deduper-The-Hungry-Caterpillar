# Deduper
### This script takes samfiles generated from single-end experiements and removes PCR duplicates.

We define PCR duplicates as reads which meet all of the following criteria:
- Have the same start position on the reference genome.
    - left most start position for forward strand alignment (after adjusting for left soft clipping).
    - right most start position for reverse strand alignment (after adjusting for total Ns, Ds, Ms, and right Ss; see cigar table (pg 8) on [sam format](https://samtools.github.io/hts-specs/SAMv1.pdf).
- Belong to the same chromosome.
- Align to the same strand.
- And have the same Unique Molecular Identifier (UMI).

## Set up
Set the environment by saving the `dedup_environment.yml` to your working directory. \
\
Then run: \
`conda env create -f dedup_environment.yml`

Then verify you have the environment via: \
`conda activate dedup`

## Input
To run the deduper script, use: `./deduper.sh <options>`
The options are as follows:

- `-h` displays help menu.
- `-f` (required) input file to be deduped. Must be an **aligned** samfile. See `-s` if not aligned.
- `-s` pass any text to this option if the input sam has not been aligned.
- `-u` input UMI file. The UMI file must contain **only** UMI sequences separated by newlines. If passed then reads not containing an UMI specified in the UMI file will be discarded.
- `-p` pass any text to this option if your samfile comes from paired end data. If you pass this option the script will happily tell you that it does not accept paired end data...yet.

## Output
The `deduper.sh` script will output:
- A sam formated file with PCR duplicates removed called `<input file name_deduped>`.
    - Note that the first encountered read will be saved and any duplicates of this first read will be removed.
- A standard out text listing:
    - The number of unique reads in each chromosome
    - The total number of reads in original sam file
    - The total number of reads in the deduped sam file
    - The total number of reads discarded due to non-valid UMI
    - The total number of PCR duplicate reads discarded
    - The percentage of removed reads relative to the original sam file

## Example usage:
In your working directory: \
```./deduper.sh -f sammy.sam -s y -u UMIs.txt```

This will first sort `sammy.sam` via [samtools](http://www.htslib.org/doc/samtools.html), and then remove PCR duplicates from `sammy.sam`, outputting the file `sammy.sam_deduped`; this out file will have all PCR duplicates reads and reads with UMIs not in `UMI.txt` removed.
