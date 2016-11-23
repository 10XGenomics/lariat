// Copyright (c) 2015 10X Genomics, Inc. All rights reserved.

/* This provides a Go <--> BWA interface */

package gobwa

// #cgo LDFLAGS: -lbwa -lz -lm -ljemalloc
// #include "bwa_bridge.h"
// #include "bwa/bwamem.h"
// #include "bwa/bwa.h"
// #include "bwa/bwt.h"
// #include <stdlib.h>
import "C"
import "unsafe"
import "log"

import "fmt"
/*
 * Holds a loaded BWA reference object.
 */
type GoBwaReference struct {
	BWTData unsafe.Pointer // Secret pointer to a *btw_t type
    contigTids map[string]int32
}

//gets contig names and lengths
func (r GoBwaReference) GetReferenceContigsInfo() ([]string, []int64) {
	typedRef := (*C.bwaidx_t)(r.BWTData)
	contigs := typedRef.bns
	numContigs := int(contigs.n_seqs)
	names := make([]string, numContigs)
	lengths := make([]int64, numContigs)
	for i := 0; i < numContigs; i++ {
		contig_ptr := (uintptr(unsafe.Pointer(contigs.anns)) + uintptr(i)*unsafe.Sizeof(*contigs.anns))
		contig := (*C.bntann1_t)(unsafe.Pointer(contig_ptr))
		names[i] = C.GoString(contig.name)
		lengths[i] = int64(int(contig.len))
	}
	return names, lengths
}

/*
 * Sets a set of BWA settings
 */
type GoBwaSettings struct {
	Settings unsafe.Pointer // Secret pointer to a *mem_opt_t type
}


func (r GoBwaReference) GetSeq(chrom string, start, end int64, reversed bool) []byte {
    typedRef := (*C.bwaidx_t)(r.BWTData)
    contigTid := r.contigTids[chrom]
    contigs := typedRef.bns
    contig_ptr := (uintptr(unsafe.Pointer(contigs.anns)) + uintptr(contigTid)*unsafe.Sizeof(*contigs.anns))
    contig := (*C.bntann1_t)(unsafe.Pointer(contig_ptr))
    offset := int64(contig.offset)
    offstart := start + offset
    offend := end + offset
    result := uintptr(unsafe.Pointer(C.bns_fetch_seq(typedRef.bns, typedRef.pac, (*C.int64_t)(&offstart), (C.int64_t)((offstart+offend)>>1), (*C.int64_t)(&offend), (*C.int)(&contigTid))))
    //len := int(end-start)
    //str := *(*[len]byte)(result)
    str := make([]byte, int(end-start))
    if reversed {
        for i := 0; i < int(offend-offstart); i++ {
            ptr := (uintptr(unsafe.Pointer(result + uintptr(i))))
            str[int(offend-offstart) - i - 1] = twoBitToSeqComp[*(*byte)(unsafe.Pointer(ptr))]
        }
    } else {
        for i := 0; i < int(offend-offstart); i++ {
            ptr := (uintptr(unsafe.Pointer(result + uintptr(i))))
            if int(*(*byte)(unsafe.Pointer(ptr))) > len(twoBitToSeq) {
                
            fmt.Println(*(*byte)(unsafe.Pointer(ptr)))
            }
            str[i] = twoBitToSeq[*(*byte)(unsafe.Pointer(ptr))]
        }
    }
    C.free(unsafe.Pointer(result));
    return str
}

var twoBitToSeq = [4]byte{'A','C','G','T'}
var twoBitToSeqComp = [4]byte{'T','G','C','A'}

/*
 * Represents a candidate alignment
 */
type EasyAlignment struct {
	Offset        int64
	Alignment_end int64
	Contig        string
	Reversed      bool
	ChainedHit    *C.mem_alnreg_t
	Score         int
	Secondary     bool
	ReadS         int
	ReadE         int
}

type Chain struct {
	Offset        int64
	Contig        string
	Reversed      bool
	chain_pointer *C.mem_chain_t
}

type Arena struct {
	Pointers []uintptr
}

func (a *Arena) Push(p uintptr) {
	a.Pointers = append(a.Pointers, p)
}

func (a *Arena) Free() {
	for i := 0; i < len(a.Pointers); i++ {
		C.free((unsafe.Pointer(a.Pointers[i])))
	}
	a.Pointers = a.Pointers[0:0]
}

func NewArena() *Arena {
	a := new(Arena)
	a.Pointers = make([]uintptr, 0, 0)
	return a
}

func GoBwaLoadReference(path string) *GoBwaReference {
	var ref *C.bwaidx_t
	ref = C.bwa_idx_load(C.CString(path), C.BWA_IDX_ALL)

	if uintptr(unsafe.Pointer(ref)) == uintptr(0) {
		log.Printf("THAT DIDN'T WORK!!! FILE PROBLEMS?")
		return nil
	}
    contigTids := map[string]int32{}
    typedRef := (*C.bwaidx_t)((unsafe.Pointer)(ref))
    contigs := typedRef.bns
    numContigs := int(contigs.n_seqs)
    for i := 0; i < numContigs; i++ {
        contig_ptr := (uintptr(unsafe.Pointer(contigs.anns)) + uintptr(i)*unsafe.Sizeof(*contigs.anns))
        contig := (*C.bntann1_t)(unsafe.Pointer(contig_ptr))
        name := C.GoString(contig.name)
        contigTids[name] = int32(i)
    }
	return &GoBwaReference{BWTData:(unsafe.Pointer)(ref), contigTids: contigTids}
}

func GoBwaAllocSettings() *GoBwaSettings {
	var s GoBwaSettings
	s.Settings = (unsafe.Pointer)(C.mem_opt_init())
	return &s
}

/*
This takes a sequence written in ASCII (AaCcTtGg) and converts it to
an array of bytes with A-->0, C-->1 G-->2, T-->3 anything else -->4
*/
func SequenceConvert(seq string) []byte {

	buffer := make([]byte, len(seq))

	for i := 0; i < len(seq); i++ {
		buffer[i] = (byte)(C.nst_nt4_table[(int)(seq[i])])
	}
	return buffer
}

/*
 * This attempts to align "seq", which is a string of ACGTacgt letters. It returns
 * an array of EasyAlignment objects.
 * TODO: Actually return that list.
 */
func GoBwaAlign(ref *GoBwaReference, settings *GoBwaSettings, seq string, arena *Arena) []EasyAlignment {

	converted_seq := SequenceConvert(seq)
	typed_ref := (*C.bwaidx_t)(ref.BWTData)

	//results := C.mem_chain((*C.mem_opt_t)(settings.Settings), typed_ref.bwt, typed_ref.bns, len(converted_seq), &(converted_seq[0]), unsafe.Pointer(uintptr(0)));

	results := C.mem_align1_core((*C.mem_opt_t)(settings.Settings),
		typed_ref.bwt,
		typed_ref.bns,
		typed_ref.pac,
		(C.int)(len(converted_seq)),
		(*C.char)(unsafe.Pointer((&(converted_seq[0])))),
		unsafe.Pointer(uintptr(0)))
	//algns := make([]*C.mem_alnreg_t, (int)(results.n))
	algns := make([]EasyAlignment, (int)(results.n))

	arena.Push(uintptr(unsafe.Pointer(results.a)))
	for i := (uintptr(0)); i < (uintptr)(results.n); i++ {
		a := ((*C.mem_alnreg_t)(unsafe.Pointer(uintptr(unsafe.Pointer(results.a)) + i*(unsafe.Sizeof(*results.a)))))
		//arena.Push(uintptr(unsafe.Pointer(a)));
		p := InterpretAlign(ref, a)
		algns[i] = p
		//log.Printf("%v", p)
		//log.Printf("%v", *((*C.mem_alnreg_t)(unsafe.Pointer(uintptr(unsafe.Pointer(results.a)) + i*(unsafe.Sizeof(*results.a))))))
	}

	return algns
	//log.Printf("%v", results);
}

func GoBwaChain(ref *GoBwaReference, settings *GoBwaSettings, seq string) []Chain {
	converted_seq := SequenceConvert(seq)
	typed_ref := (*C.bwaidx_t)(ref.BWTData)

	results := C.mem_chain((*C.mem_opt_t)(settings.Settings),
		typed_ref.bwt,
		typed_ref.bns,
		(C.int)(len(converted_seq)),
		(*C.uint8_t)(&(converted_seq[0])),
		unsafe.Pointer(uintptr(0)))
	chns := make([]Chain, (int)(results.n))

	for i := (uintptr(0)); i < (uintptr)(results.n); i++ {
		a := ((*C.mem_chain_t)(unsafe.Pointer(uintptr(unsafe.Pointer(results.a)) + i*(unsafe.Sizeof(*results.a)))))
		p := InterpretChain(ref, a)
		chns[i] = p
		//log.Printf("%v", p)
	}
	return chns
}

func GoBwaMemMateSW(ref *GoBwaReference, settings *GoBwaSettings, read1 *[]byte, read2 *[]byte, arena *Arena, score_delta int) ([]EasyAlignment, []EasyAlignment) {

	typed_ref := (*C.bwaidx_t)(ref.BWTData)
	var Pes [4]C.mem_pestat_t
	Pes[0].failed = 1
	Pes[1].low = -35
	Pes[1].high = 500
	Pes[1].failed = 0
	Pes[1].avg = 200.0
	Pes[1].std = 100.0
	Pes[2].failed = 1
	Pes[3].failed = 1
	converted_seq_read1 := SequenceConvert(string(*read1))
	converted_seq_read2 := SequenceConvert(string(*read2))
	// get mapping for each read
    read1_results := C.mem_alnreg_v{}
    read2_results := C.mem_alnreg_v{}
    if len(converted_seq_read1) > 0 {
	    read1_results = C.mem_align1_core((*C.mem_opt_t)(settings.Settings),
		    typed_ref.bwt,
		    typed_ref.bns,
		    typed_ref.pac,
		    (C.int)(len(converted_seq_read1)),
		    (*C.char)(unsafe.Pointer((&(converted_seq_read1[0])))),
		    unsafe.Pointer(uintptr(0)))
    }
    if len(converted_seq_read2) > 0 {
	    read2_results = C.mem_align1_core((*C.mem_opt_t)(settings.Settings),
		    typed_ref.bwt,
		    typed_ref.bns,
		    typed_ref.pac,
		    (C.int)(len(converted_seq_read2)),
		    (*C.char)(unsafe.Pointer((&(converted_seq_read2[0])))),
	    	unsafe.Pointer(uintptr(0)))
    }

	algns_read1 := make([]EasyAlignment, (int)(read1_results.n))
	algns_read2 := make([]EasyAlignment, (int)(read2_results.n))
	best_read1_score := 0
	best_read2_score := 0
	//get some interpretation of these alignments and note the best score for each
	for i := (uintptr(0)); i < (uintptr)(read1_results.n); i++ {
		a := ((*C.mem_alnreg_t)(unsafe.Pointer(uintptr(unsafe.Pointer(read1_results.a)) + i*(unsafe.Sizeof(*read1_results.a)))))
		p := InterpretAlign(ref, a)
		algns_read1[i] = p
		if p.Score > best_read1_score {
			best_read1_score = p.Score
		}
	}

	for i := (uintptr(0)); i < (uintptr)(read2_results.n); i++ {
		a := ((*C.mem_alnreg_t)(unsafe.Pointer(uintptr(unsafe.Pointer(read2_results.a)) + i*(unsafe.Sizeof(*read2_results.a)))))
		p := InterpretAlign(ref, a)
		algns_read2[i] = p
		if p.Score > best_read2_score {
			best_read2_score = p.Score
		}
	}

	//rescue alignments for read1 by looping through read2's hits and doing mem_matesw
	num := 0
	for i := 0; i < len(algns_read2) && num < 50 && len(converted_seq_read1) > 0; i++ {
		if algns_read2[i].Score >= best_read2_score-score_delta {
			// attempt to rescue the read1 alignment here
			num++
			C.mem_matesw((*C.mem_opt_t)(settings.Settings),
		    	typed_ref.bns,
		    	typed_ref.pac,
		    	&Pes[0],
				algns_read2[i].ChainedHit,
			    (C.int)(len(*read1)),
		    	(*C.uint8_t)(&(converted_seq_read1[0])),
		    	&read1_results)
		}
	}

	algns_read1 = make([]EasyAlignment, (int)(read1_results.n))
	for i := (uintptr(0)); i < (uintptr)(read1_results.n); i++ {
		a := ((*C.mem_alnreg_t)(unsafe.Pointer(uintptr(unsafe.Pointer(read1_results.a)) + i*(unsafe.Sizeof(*read1_results.a)))))
		p := InterpretAlign(ref, a)
		algns_read1[i] = p
	}

	//rescue alignments for read2 by looping through read1's hits and doing mem_matesw
	num = 0
	for i := 0; i < len(algns_read1) && num < 50 && len(converted_seq_read2) > 0; i++ {
		if algns_read1[i].Score >= best_read1_score-score_delta {
			// attempt to rescue the read2 alignment here
			num++
			C.mem_matesw((*C.mem_opt_t)(settings.Settings),
				typed_ref.bns,
				typed_ref.pac,
				&Pes[0],
				algns_read1[i].ChainedHit,
				(C.int)(len(*read2)),
				(*C.uint8_t)(&(converted_seq_read2[0])),
				&read2_results)
		}
	}

	arena.Push(uintptr(unsafe.Pointer(read1_results.a)))
	arena.Push(uintptr(unsafe.Pointer(read2_results.a)))

	algns_read2 = make([]EasyAlignment, (int)(read2_results.n))

	for i := (uintptr(0)); i < (uintptr)(read2_results.n); i++ {
		a := ((*C.mem_alnreg_t)(unsafe.Pointer(uintptr(unsafe.Pointer(read2_results.a)) + i*(unsafe.Sizeof(*read2_results.a)))))
		p := InterpretAlign(ref, a)
		algns_read2[i] = p
	}
	return algns_read1, algns_read2
}

func InterpretAlign(ref *GoBwaReference, caln *C.mem_alnreg_t) EasyAlignment {
	var res EasyAlignment

	typed_ref := (*C.bwaidx_t)(ref.BWTData)
	contigs := typed_ref.bns

	var contig_id = uintptr(caln.rid)

	contig_ptr := (uintptr(unsafe.Pointer(contigs.anns)) + contig_id*unsafe.Sizeof(*contigs.anns))

	contig := (*C.bntann1_t)(unsafe.Pointer(contig_ptr))
	// log.Printf("CONTIG PTR: %v ", caln)
	// log.Printf("what is it that we want %p %v %v ",caln, *caln, contig_id)
	if caln.rb < contigs.l_pac {
		res.Offset = int64(caln.rb - contig.offset)
		res.Reversed = false
	} else {
		res.Offset = int64(contigs.l_pac*2 - 1 - caln.rb - contig.offset)
		res.Reversed = true
	}
	if caln.re < contigs.l_pac {
		res.Alignment_end = int64(caln.re) - int64(contig.offset)
	} else {
		res.Alignment_end = int64(contigs.l_pac*2 - 1 - caln.re - contig.offset)
	}
	res.Contig = C.GoString(contig.name)
	res.Secondary = int(caln.secondary) >= 0 || int(caln.secondary_all) > 0
	res.ChainedHit = caln
	res.Score = int(caln.score)
	res.ReadS = int(caln.qb)
	res.ReadE = int(caln.qe)
	return res
}

func InterpretChain(ref *GoBwaReference, chn *C.mem_chain_t) Chain {
	var res Chain

	typed_ref := (*C.bwaidx_t)(ref.BWTData)
	contigs := typed_ref.bns

	var contig_id = uintptr(chn.rid)

	contig_ptr := (uintptr(unsafe.Pointer(contigs.anns)) + contig_id*unsafe.Sizeof(*contigs.anns))

	contig := (*C.bntann1_t)(unsafe.Pointer(contig_ptr))
	log.Printf("what is it that we want %v %v %v ", chn, *chn)
	log.Printf("CONTIG PTR: %v %v %v %v", contig, contig.offset, contig.name, contig.anno)

	if chn.pos < contigs.l_pac {
		res.Offset = int64(chn.pos - contig.offset)
		res.Reversed = false
	} else {
		res.Offset = int64(contigs.l_pac*2 - 1 - chn.pos - contig.offset)
		res.Reversed = true
	}

	res.Contig = C.GoString(contig.name)
	res.chain_pointer = chn
	return res
}

func GoBwaSmithWaterman(ref *GoBwaReference, settings *GoBwaSettings, seq string, alignment unsafe.Pointer, arena *Arena) SingleReadAlignment {
	converted_seq := SequenceConvert(seq)
	typed_alignment := (*C.mem_alnreg_t)(alignment)
	typed_ref := (*C.bwaidx_t)(ref.BWTData)
	results := C.mem_reg2aln((*C.mem_opt_t)(settings.Settings),
		typed_ref.bns,
		typed_ref.pac,
		(C.int)(len(converted_seq)),
		(*C.char)(unsafe.Pointer((&(converted_seq[0])))),
		typed_alignment)
	// arena.Push(uintptr(unsafe.Pointer(results)));
	arena.Push(uintptr(unsafe.Pointer(results.cigar)))
	arena.Push(uintptr(unsafe.Pointer(results.XA)))

	return InterpretSingleReadAlignment(ref, &results)
}

type SingleReadAlignment struct {
	Pos                 int64
	Chrom               string
	Flag                int //what is this?
	Reversed            bool
	Alt                 int
	Mapq                int
	EditDistance        int
	Cigar               []uint32
	AlternativeMappings string
	Score               int
	Sub                 int //what is this?
	AltSC               int //what is this?
	ReadS               int //What part of the read is covered by this alignment
	ReadE               int
	raw                 *C.mem_aln_t
}

func EnumerateContigs(ref *GoBwaReference, callback func(name string, length int)) {
	typed_ref := (*C.bwaidx_t)(ref.BWTData)

	contigs := (typed_ref.bns)

	for i := uintptr(0); i < uintptr(contigs.n_seqs); i++ {
		cptr := (uintptr(unsafe.Pointer(contigs.anns)) + i*unsafe.Sizeof(*contigs.anns))
		c := (*C.bntann1_t)(unsafe.Pointer(cptr))

		callback(C.GoString(c.name), int(c.len))
	}

}

func InterpretSingleReadAlignment(ref *GoBwaReference, alignment *C.mem_aln_t) SingleReadAlignment {
	var result SingleReadAlignment

	fixed_flags := ((*C.GO_mem_aln_t)((unsafe.Pointer)(alignment))).flag2
	typed_ref := (*C.bwaidx_t)(ref.BWTData)
	contigs := typed_ref.bns
	var contig_id = uintptr(alignment.rid)
	contig_ptr := (uintptr(unsafe.Pointer(contigs.anns)) + contig_id*unsafe.Sizeof(*contigs.anns))
	contig := (*C.bntann1_t)(unsafe.Pointer(contig_ptr))
	result.Pos = int64(alignment.pos)
	result.Chrom = C.GoString(contig.name)
	cigar_ops := uintptr(alignment.n_cigar)
	raw_cigar := make([]uint32, cigar_ops)
	//fmt.Println(result.Chrom)
	//fmt.Println(result.Pos)
	result.Cigar = make([]uint32, cigar_ops*2)
	for i := uintptr(0); i < cigar_ops; i++ {
		raw_cigar[i] = *(*uint32)(unsafe.Pointer(uintptr(unsafe.Pointer(alignment.cigar)) + i*unsafe.Sizeof(C.uint32_t(0))))
		//fmt.Println(raw_cigar[i]>>4)
		//fmt.Println(string([]byte{'M','I','D','S','H'}[raw_cigar[i]&0xf]))
		result.Cigar[i*2] = raw_cigar[i] & 0xf
		result.Cigar[i*2+1] = raw_cigar[i] >> 4
	}
	//log.Printf("Fixed flags %x", int(fixed_flags))
	result.Alt = int(fixed_flags & 0x2)
	//fmt.Println("is alt ",result.Alt)
	result.Mapq = int(fixed_flags&0x2C) >> 2
	result.Reversed = (uint32(0) != uint32(uint32(fixed_flags)&uint32(0x1)))
	//fmt.Println("mapq ",result.Mapq)
	result.Score = int(alignment.score)
	//fmt.Println("score ",result.Score)
	result.Sub = int(alignment.sub)
	//fmt.Println("sub ",result.Sub)
	result.EditDistance = int(uint(fixed_flags) >> 10)
	// fmt.Println("NM ",int(alignment.flag&0x3FFFFF))
	result.AltSC = int(alignment.alt_sc)
	//fmt.Println("alt sc ",result.AltSC)
	result.raw = alignment
	return result
}
