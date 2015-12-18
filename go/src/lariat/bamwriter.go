// Copyright (c) 2015 10X Genomics, Inc. All rights reserved.

package main

import (
	"code.google.com/p/biogo.bam"
	"fmt"
	"gobwa"
	"log"
	"os"
	"strconv"
	"strings"
	"sync"
	"time"
)

const (
	SAM_CIGAR_MATCH     = 0
	SAM_CIGAR_INSERT    = 1
	SAM_CIGAR_DEL       = 2
	SAM_CIGAR_SKIP      = 3
	SAM_CIGAR_SOFT_CLIP = 4
	SAM_CIGAR_HARD_CLIP = 5
)

type BAMWriters struct {
	BarcodeSortedBam     *BAMWriter
	PositionBucketedBams map[string][]*BAMWriter
	positionChunkSize    int
	debugTags            bool
	channel              chan *Data
	/* This mutex is Rlocked by the worker thread. When we close, we wait
	 * for the mutex to be unlocked to ensure data is flushed before continueing
	 */
	done sync.RWMutex
}

type BAMWriter struct {
	Writer  *bam.Writer
	Contigs map[string]*bam.Reference
	Record  bam.Record
}

func CreateBAM(ref *gobwa.GoBwaReference, path, read_group, sample_id string) (*BAMWriter, error) {
	bw := &BAMWriter{}
	bw.Contigs = make(map[string]*bam.Reference)

	references := make([]*bam.Reference, 0, 0)

	gobwa.EnumerateContigs(ref, func(name string, length int) {
		r, err := bam.NewReference(name, name, "human", length, nil, nil)
		if err != nil {
			panic(err)
		}
		references = append(references, r)
		bw.Contigs[name] = r
	})

	// NewReadGroup(name, center, desc, lib, prog, plat, unit, sample string, date time.Time, size int, flow, key []byte)
	rg, err := bam.NewReadGroup(
		read_group,
		"",
		"",
		"",
		"",
		"",
		"",
		sample_id,
		time.Now(),
		0,
		nil,
		nil)

	if err != nil {
		panic(err)
	}

	h, err := bam.NewHeader(nil, references)

	if err != nil {
		panic(err)
	}

	h.AddReadGroup(rg)

	// Add a program line for lariat
	prog := bam.NewProgram("lariat", "lariat", strings.Join(os.Args, " "), "", __VERSION__)
	h.AddProgram(prog)

	file, err := os.Create(path)

	if err != nil {
		return nil, err
	}

	w, err := bam.NewWriter(file, h, 2)

	if err != nil {
		panic(err)
	}
	bw.Writer = w
	return bw, nil
}

func (b BAMWriters) getPositionBucketedBamForAlignment(aln *Alignment, unmapped bool) *BAMWriter {
	if unmapped {
		return b.PositionBucketedBams["unmapped"][0]
	}
	return b.PositionBucketedBams[aln.contig][aln.pos/int64(b.positionChunkSize)]
}

func CreateBAMs(ref *gobwa.GoBwaReference, basePath, read_group, sample_id string, positionChunkSize int, debugTags bool) (*BAMWriters, error) {
	barcodeSortedBam, err := CreateBAM(ref, basePath+"/bc_sorted_bam.bam", read_group, sample_id)
	if err != nil {
		return nil, err
	}
	contigNames, contigLengths := ref.GetReferenceContigsInfo()
	PositionBucketedBams := make(map[string][]*BAMWriter, len(contigNames))

	for index, contigName := range contigNames {
		offset := int64(0)
		PositionBucketedBams[contigName] = make([]*BAMWriter, contigLengths[index]/int64(positionChunkSize)+int64(1))
		for chunkIndex := 0; offset < contigLengths[index]; chunkIndex, offset = chunkIndex+1, offset+int64(positionChunkSize) {
			offsetString := strconv.FormatInt(offset, 10)
			newOffsetString := strconv.FormatInt(offset, 10)
			for digit := 0; digit < 10; digit++ {
				if digit < len(offsetString) {
					continue
				}
				newOffsetString = "0" + newOffsetString
			}
			leadingZeros := ""
			for digit := 0; digit < 15-len(strconv.Itoa(index)); digit++ {
				leadingZeros += "0"
			}
			PositionBucketedBams[contigName][chunkIndex], err = CreateBAM(ref, basePath+"/"+leadingZeros+strconv.Itoa(index)+"-"+contigName+"_"+newOffsetString+"_pos_bucketed.bam", read_group, sample_id)
			if err != nil {
				return nil, err
			}
		}
	}
	unmappedBam, err := CreateBAM(ref, basePath+"/"+"ZZZ_unmapped_pos_bucketed.bam", read_group, sample_id)
	if err != nil {
		return nil, err
	}
	PositionBucketedBams["unmapped"] = []*BAMWriter{unmappedBam}
	toReturn := BAMWriters{BarcodeSortedBam: barcodeSortedBam, PositionBucketedBams: PositionBucketedBams, positionChunkSize: positionChunkSize, debugTags: debugTags}
	toReturn.channel = make(chan *Data, 8)
	go BamThread(&toReturn)
	return &toReturn, nil
}

func auxify_string(name []byte, data []byte) []byte {
	vec := make([]byte, len(data)+3)

	vec[0] = name[0]
	vec[1] = name[1]
	vec[2] = byte('Z')
	for i := 0; i < len(data); i++ {
		vec[3+i] = data[i]
	}
	return vec
}

func FixCigar(in []uint32) bam.Cigar {

	count := (len(in) / 2)
	cigar := make(bam.Cigar, count)

	for i := (0); i < count; i++ {
		cigar[i] = bam.NewCigarOp(bam.CigarOpType(in[i*2]), int(in[i*2+1]))
	}

	return cigar
}

func fixQual(in []byte) []byte {
	output := make([]byte, len(in))
	for i := 0; i < len(in); i++ {
		output[i] = in[i] - 33
	}

	return output
}

var cigartable = [5]uint32{
	0: 0,
	1: 1,
	2: 2,
	3: 4,
	4: 5,
}

func fixCigar(in []uint32) []uint32 {
	var out = make([]uint32, len(in))

	for i := 0; i < len(in)/2; i++ {
		idx := i * 2
		if int(in[idx]) > len(cigartable) {
			log.Printf("BAMOP: %v", in[idx])
			panic("ILLEGAL CIGAR OP")
		}
		out[idx] = cigartable[int(in[idx])]
		out[idx+1] = in[idx+1]
	}
	return out
}

func (b *BAMWriters) AppendBams(aln *Alignment, primary *Alignment, debugTags bool, attach_bx bool) {
	b.BarcodeSortedBam.AppendBam(aln, primary, debugTags, attach_bx)
	b.getPositionBucketedBamForAlignment(aln, aln.IsUnmapped()).AppendBam(aln, primary, debugTags, attach_bx)
}

func (b *BAMWriter) AppendBam(aln *Alignment, primary *Alignment, debugTags bool, attach_bx bool) {
	ref := b.Contigs[aln.contig]
	var flags int32

	if !aln.is_proper && aln.score-17 < 19 {
		aln.pos = -1
		aln.mapq = 0
	}
	if aln.mate_id >= 0 {
		flags |= 1
		if aln.is_proper {
			flags |= 0x2
		}
		if primary.mate_alignment.pos == -1 || (!primary.is_proper && primary.mate_alignment.score-17 < 19) {
			flags |= 0x8
			b.Record.MatePos = -1
			b.Record.MateRef = nil
		}
		if primary.mate_alignment.reversed {
			flags |= 0x20
		}
		if aln.read1 {
			flags |= 0x40
		} else {
			flags |= 0x80
		}
		if aln.duplicate {
			flags |= 0x400
		}

		b.Record.MateRef = b.Contigs[primary.mate_alignment.contig]
		b.Record.MatePos = int(primary.mate_alignment.pos)
		if primary.mate_alignment.pos == -1 {
			b.Record.MateRef = nil
		} else if aln == primary {
			if aln.contig == aln.mate_alignment.contig {
				if aln.reversed {
					b.Record.TempLen = -int(aln.aend - aln.mate_alignment.pos)
				} else {
					b.Record.TempLen = int(aln.mate_alignment.aend - aln.pos)
				}
			}

		}

	} else {
		b.Record.MatePos = -1
		b.Record.MateRef = nil
	}

	if aln != primary {
		flags |= 256
	}

	b.Record.Ref = ref

	b.Record.MapQ = byte(aln.mapq)
	if aln.pos == -1 {
		flags |= 0x4
		b.Record.MapQ = byte(0)
		b.Record.Ref = nil
	}
	if aln.reversed {
		flags |= 0x10
	}
	b.Record.Name = strings.TrimRight(*(aln.read_name), "\n")

	b.Record.Flags = bam.Flags(flags)

	seq := *aln.read_seq
	pos := int(aln.pos)
	cigar := fixCigar(aln.cigar)
	qual := *aln.read_qual

	if aln.reversed {
		seq = reverseComp(seq)
		qual = reverseQual(qual)
	}

	if primary != aln {
		var deltapos int
		seq, qual, cigar = HardClip(seq, qual, cigar, aln.reversed)
		if pos > 0 {
			pos += deltapos
		}
	}

	b.Record.Pos = pos
	b.Record.Cigar = FixCigar(cigar)
	b.Record.Seq = bam.NewNybbleSeq(seq)
	b.Record.Qual = fixQual(qual)

	barcode := strings.Split(string(*aln.raw_barcode), "-")
	bwa_pick := "0"
	if aln.bwa_pick {
		bwa_pick = "1"
	}
    
	aux := []bam.Aux{}
	bw := auxify_string([]byte("BWA"), []byte(bwa_pick))
    aux = append(aux, bam.Aux(bw))
    if len(barcode[0]) > 1 {
	    rx := auxify_string([]byte("RX"), []byte(barcode[0]))
	    qx := auxify_string([]byte("QX"), *aln.barcode_qual)
        aux = append(aux, bam.Aux(rx))
        aux = append(aux, bam.Aux(qx))
    }
	bc := auxify_string([]byte("BC"), *aln.sample_index)
	qt := auxify_string([]byte("QT"), *aln.sample_index_qual)
	rg := auxify_string([]byte("RG"), []byte(*aln.read_group))
	as := auxify_string([]byte("AS"), []byte(strconv.FormatInt(int64(aln.score), 10)))
	xs := auxify_string([]byte("XS"), []byte(strconv.FormatInt(int64(0), 10)))
    if len(*aln.sample_index) > 1 {
        aux = append(aux, bam.Aux(bc))
        aux = append(aux, bam.Aux(qt))
    }
    aux = append(aux, bam.Aux(as))
    aux = append(aux, bam.Aux(rg))
    aux = append(aux, bam.Aux(xs))
	if aln.mapq_data != nil && aln.mapq_data.second_best != nil {
		xs = auxify_string([]byte("XS"), []byte(strconv.FormatInt(int64(aln.mapq_data.second_best.score), 10)))
	}


	if debugTags && aln.mapq_data != nil {
		cp := auxify_string([]byte("CP"), []byte(strconv.FormatInt(int64(aln.mapq_data.copies), 10)))
		cm := auxify_string([]byte("CM"), []byte(strconv.FormatInt(int64(aln.mapq_data.copies_in_active_molecules), 10)))
		cu := auxify_string([]byte("CU"), []byte(strconv.FormatInt(int64(aln.mapq_data.unique_molecules_active), 10)))
		cs := auxify_string([]byte("CS"), []byte(strconv.FormatInt(int64(aln.mapq_data.copies_outside_active_molecules), 10)))
		rd := auxify_string([]byte("RD"), []byte(strconv.FormatInt(int64(aln.mapq_data.reads_in_molecule), 10)))
		pp := auxify_string([]byte("PP"), []byte(strconv.FormatBool(aln.is_proper)))
		aa := auxify_string([]byte("AA"), []byte(aln.mapq_data.active_alignments_in_molecules))
		mc := auxify_string([]byte("MC"), []byte(strconv.FormatFloat(float64(aln.molecule_confidence), 'f', 6, 64)))
		ps := auxify_string([]byte("PS"), []byte(strconv.FormatInt(int64(primary.mate_alignment.score), 10)))
		pl := auxify_string([]byte("PL"), []byte(strconv.FormatFloat(float64(primary.mate_alignment.log_alignment_probability), 'f', 6, 64)))
		ac := auxify_string([]byte("AC"), []byte("Match:"+strconv.FormatInt(int64(aln.matches), 10)+":Mismatches:"+strconv.FormatInt(int64(aln.mismatches), 10)+":Indels:"+strconv.FormatInt(int64(aln.indels), 10)+":soft_clipped:"+strconv.FormatBool(aln.soft_clipped)))
		pc := auxify_string([]byte("PC"), []byte("Match:"+strconv.FormatInt(int64(primary.mate_alignment.matches), 10)+":Mismatches:"+strconv.FormatInt(int64(primary.mate_alignment.mismatches), 10)+":Indels:"+strconv.FormatInt(int64(primary.mate_alignment.indels), 10)+":soft_clipped:"+strconv.FormatBool(primary.mate_alignment.soft_clipped)))
		if aln.mapq_data.second_best != nil {
			second_best_log_probability := auxify_string([]byte("XL"), []byte(strconv.FormatFloat(aln.mapq_data.second_best.log_alignment_probability, 'f', 6, 64)))
			second_best_proper_pair := auxify_string([]byte("XP"), []byte(strconv.FormatBool(aln.mapq_data.second_best_proper_pair)))
			second_best_molecule_reads := auxify_string([]byte("XR"), []byte(strconv.FormatInt(int64(aln.mapq_data.second_best_molecule_reads), 10)))
			second_best_molecule_confidence := auxify_string([]byte("XC"), []byte(strconv.FormatFloat(aln.mapq_data.second_best_molecule_confidence, 'f', 6, 64)))
			if aln.mapq_data.second_best.mate_alignment != nil {
				xm := auxify_string([]byte("XM"), []byte(strconv.FormatFloat(float64(aln.mapq_data.second_best.mate_alignment.log_alignment_probability), 'f', 6, 64)))
				aux = append(aux, bam.Aux(xm))
				xz := auxify_string([]byte("XZ"), []byte("Match:"+strconv.FormatInt(int64(aln.mapq_data.second_best.mate_alignment.matches), 10)+":Mismatches:"+strconv.FormatInt(int64(aln.mapq_data.second_best.mate_alignment.mismatches), 10)+":Indels:"+strconv.FormatInt(int64(aln.mapq_data.second_best.mate_alignment.indels), 10)+":soft_clipped:"+strconv.FormatBool(aln.mapq_data.second_best.mate_alignment.soft_clipped)))
				aux = append(aux, bam.Aux(xz))
			}
			xx := auxify_string([]byte("XX"), []byte("Match:"+strconv.FormatInt(int64(aln.mapq_data.second_best.matches), 10)+":Mismatches:"+strconv.FormatInt(int64(aln.mapq_data.second_best.mismatches), 10)+":Indels:"+strconv.FormatInt(int64(aln.mapq_data.second_best.indels), 10)+":soft_clipped:"+strconv.FormatBool(aln.mapq_data.second_best.soft_clipped)))
			aux = append(aux, bam.Aux(xx))
			aux = append(aux, bam.Aux(second_best_log_probability))
			aux = append(aux, bam.Aux(second_best_proper_pair))
			aux = append(aux, bam.Aux(second_best_molecule_reads))
			aux = append(aux, bam.Aux(second_best_molecule_confidence))
		}
		aux = append(aux, bam.Aux(aa))
		aux = append(aux, bam.Aux(cp))
		aux = append(aux, bam.Aux(cm))
		aux = append(aux, bam.Aux(cu))
		aux = append(aux, bam.Aux(cs))
		aux = append(aux, bam.Aux(rd))
		aux = append(aux, bam.Aux(mc))
		aux = append(aux, bam.Aux(pp))
		aux = append(aux, bam.Aux(ps))
		aux = append(aux, bam.Aux(pl))
		aux = append(aux, bam.Aux(ac))
		aux = append(aux, bam.Aux(pc))
	}
	if len(barcode) > 1 && attach_bx {
		bx := auxify_string([]byte("BX"), *aln.raw_barcode) //this is confusing because it isnt the raw barcode, should rename
		aux = append(aux, bam.Aux(bx))
		if aln.active_molecule {
			md := auxify_string([]byte("DM"), []byte(strconv.FormatFloat(aln.molecule_difference, 'f', 6, 64)))
			aux = append(aux, bam.Aux(md))
		}
	}
	b.Record.AuxTags = aux
	b.Writer.Write(&b.Record)
}

func (b *BAMWriters) Close() {
	b.channel <- nil
	b.done.Lock()
	b.done.Unlock()
}

var complement = [256]byte{
	'A': 'T',
	'a': 'T',
	'C': 'G',
	'c': 'G',
	'G': 'C',
	'g': 'C',
	'T': 'A',
	't': 'A',
	'N': 'N',
	'n': 'N',
}

func reverseComp(seq []byte) []byte {
	toReturn := make([]byte, len(seq))
	for i := 0; i < len(seq); i++ {
		toReturn[i] = complement[seq[len(seq)-i-1]]
	}
	return toReturn
}

func reverseCigar(cig []uint32) []uint32 {
	toReturn := make([]uint32, len(cig))
	for i := 0; i < len(cig); i += 2 {
		toReturn[i+1] = cig[len(cig)-i-1]
		toReturn[i] = cig[len(cig)-i-2]
	}
	return toReturn
}

func reverseQual(qual []byte) []byte {
	toReturn := make([]byte, len(qual))
	for i := 0; i < len(qual); i++ {
		toReturn[i] = qual[len(qual)-i-1]
	}
	return toReturn
}

func DumpToBams(alignments *Data, b *BAMWriters) {
	b.channel <- alignments
}

func BamThread(b *BAMWriters) {
	b.done.RLock()

	for alignments := <-b.channel; alignments != nil; alignments = <-b.channel {
		DoDumpToBam(alignments.alignments, b, b.debugTags, alignments.attach_bx)
		ReturnBuffer(alignments.reads)
	}
	b.BarcodeSortedBam.Writer.Close()
	for _, bams := range b.PositionBucketedBams {
		for _, bam := range bams {
			bam.Writer.Close()
		}
	}
	b.done.RUnlock()
}

func DoDumpToBam(alignments [][]*Alignment, b *BAMWriters, debugTags bool, attach_bx bool) {
	reads := 0
	for _, alignmentArray := range alignments {
		if alignmentArray == nil || len(alignmentArray) == 0 {
			panic(fmt.Sprintf("not all read_ids are spoken for"))
		}
		read_output := false
		if alignmentArray != nil {
			for _, alignment := range alignmentArray {
				if alignment.active {
					b.AppendBams(alignment, alignment, debugTags, attach_bx)
					if alignment.secondary != nil {
						b.AppendBams(alignment.secondary, alignment, debugTags, attach_bx)
					}
					reads++
					read_output = true
				}
			}
		}
		if !read_output {
			panic(fmt.Sprintf("read_id has no active alignment but more than one alignment"))
		}
	}
}

/*
 Convert from "soft" clipping to "hard" clipping. Truncate the sequence and quality
 and convert "S" to "H" in teh cigar string.
*/
func HardClip(seq []byte, qual []byte, cigar []uint32, reversed bool) ([]byte, []byte, []uint32) {
	var start, end int

	start = 0
	end = len(seq)

	newcigar := make([]uint32, len(cigar))
	copy(newcigar, cigar)
	if len(newcigar) >= 2 {
		if newcigar[0] == SAM_CIGAR_SOFT_CLIP {
			start = int(newcigar[1])
			newcigar[0] = SAM_CIGAR_HARD_CLIP
		}
	}
	if len(newcigar) >= 4 {
		p := len(newcigar) - 2
		if newcigar[p] == SAM_CIGAR_SOFT_CLIP {
			end -= int(newcigar[p+1])
			newcigar[p] = SAM_CIGAR_HARD_CLIP
		}
	}

	newseq := (seq)[start:end]
	newqual := (qual)[start:end]
	return newseq, newqual, newcigar
}
