# Lariat: Linked-Read Alignment Tool

Lariat is an aligner for barcoded linked reads, produced by the 10X Genomics GemCode™ platform. All the linked reads for a single barcode are aligned simultaneously, with the prior knowledge that the reads arise from a small number of long (10kb - 200kb) molecules. This approach allows reads to be mapped in repetitive regions of the genome.

We follow the RFA method developed in the Batzoglou lab by Bishara, Liu et. al. [Genome Res. 2015. 25:1570-1580](http://genome.cshlp.org/content/25/10/1570) 

Lariat generates candidate alignments by calling the BWA C API, then performs the RFA inference to select the final mapping position and MAPQ.

## Usage Notes: 
* Lariat currently requires Go 1.3
* Lariat currently requires a non standard format for input reads. We recommend using the Lariat build bundled with the 10X Genomics Longranger software (http://software.10xgenomics.com/)

Please contact us if you're interested in using Lariat independently of the Longranger pipeline.

## Input File Format

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
Lariat is distributed under the [GPLv3](http://www.gnu.org/licenses/gpl-3.0.en.html). Lariat links to [BWA](https://github.com/lh3/bwa) at the object level. Lariat include the BWA source code via git submodule. BWA is licensed also licensed under the GPLv3. 
