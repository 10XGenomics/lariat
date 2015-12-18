# Lariat: Linked-Read Alignment Tool

Lariat is an aligner for barcoded linked reads, produced by the 10X Genomics GemCode™ platform. All the linked reads for a single barcode are aligned simultaneously, with the prior knowledge that the reads arise from a small number of long (10kb - 200kb) molecules. This approach allows reads to be mapped in repetitive regions of the genome.

We follow the RFA method developed in the Batzoglou lab by Bishara, Liu et. al. [Genome Res. 2015. 25:1570-1580](http://genome.cshlp.org/content/25/10/1570) 

Lariat generates candidate alignments by calling the BWA C API, then performs the RFA inference to select the final mapping position and MAPQ.

## License
Lariat is distributed under the [GPLv3](http://www.gnu.org/licenses/gpl-3.0.en.html). Lariat links to [BWA](https://github.com/lh3/bwa) at the object level. Lariat include the BWA source code via git submodule. BWA is licensed also licensed under the GPLv3. 
