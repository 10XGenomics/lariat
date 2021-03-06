# Lariat: Linked-Read Alignment Tool

Lariat is an aligner for barcoded linked reads, produced by the 10X Genomics GemCode™ platform. All the linked reads for a single barcode are aligned simultaneously, with the prior knowledge that the reads arise from a small number of long (10kb - 200kb) molecules. This approach allows reads to be mapped in repetitive regions of the genome.

Lariat is based on the original RFA method developed by Alex Bishara, Yuling Liu et al in Serafim Batzoglou’s lab at Stanford: [Genome Res. 2015. 25:1570-1580](http://genome.cshlp.org/content/25/10/1570).  In addition to developing the original model for RFA, Alex Bishara and Yuling Liu both contributed substantially to the Lariat implementation maintained in this repository.

Lariat generates candidate alignments by calling the BWA C API, then performs the RFA inference to select the final mapping position and MAPQ.

## Usage Notes: 

*NOTE*: If you just want to get Lariat-aligned BAM files from Chromium Linked-Read data, you can run the ALIGN pipeline in [Long Ranger 2.2](https://support.10xgenomics.com/genome-exome/software/downloads/latest). It runs the FASTQ processing and alignment steps only.


* Lariat currently is tested with Go version 1.9.2.
* Lariat currently requires a non standard format for input reads. We recommend using the Lariat build bundled with the 10X Genomics Long Ranger software (http://software.10xgenomics.com/)

Please contact us if you're interested in using Lariat independently of the Long Ranger pipeline.

## Build notes:
In the lariat directory, run `git submodule --init --recursive` to ensure you've checked out the BWA submodule.

Make sure you have a working Go installation (version >= 1.9.2). `go version` should return something like "go version go1.9.2 linux/amd64"

From the root of the repo:
```
cd go
make           # Build lariat
bin/lariat -h  # Show cmd-line flags
```

For experimental purposes you can replace the lariat binary in a Long Ranger build with bin/lariat.


## Input File Format

The SORT_FASTQS stage in Long Ranger creates specially formatted, barcode sorted input for lariat.  We recommend using those input files to experiment with changes to lariat.
Lariat requires input data in a non-standard FASTQ-like format. Each read-pair is formatted as a record of 9 consecutive lines containing:
* read header
* read1 sequence
* read1 quals
* read2 sequence
* read2 quals
* 10X barcode string
* 10X barcode quals
* sample index sequence
* sample index quals

Read pairs must be sorted by the 10X barcode string. The 10X barcode string is of the form 'ACGTACGTACGTAC-1'. 

## License
Lariat is distributed under the MIT license. Lariat links to [BWA](https://github.com/lh3/bwa) at the object level. Lariat include the BWA source code via git submodule. Lariat links to the Apache2 branch of the BWA repo, which is licensed under the Apache2 license.
