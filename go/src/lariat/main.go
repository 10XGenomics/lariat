// Copyright (c) 2015 10X Genomics, Inc. All rights reserved.

package main

import (
	"C"
	"bufio"
	"crypto/md5"
	"encoding/binary"
	"fastqreader"
	"flag"
	"fmt"
	. "gobwa"
	"math"
	"math/rand"
	"optimizer"
	"os"
	"runtime"
	"sort"
	"strconv"
	"strings"
	"sync"
	"sync/atomic"
	"syscall"
	"unsafe"
)

// build version -- get set statically by a linker flag
var __VERSION__ string

/*Command line arguments*/
var reads = flag.String("reads", "", "fastq.gz input file containing reads (required)")
var improper_pair_penalty = flag.Float64("improper_pair_penalty", -17, "Penalty for improper pair")
var molecule_penalty = flag.Float64("molecule_penalty", -25, "Penalty for opening each molecule")
var log10_inactive_molecule = flag.Float64("log10_inactive_molecule", -3.7, "alignment log10 probability of being outside an active molecule when an alignment inside an active molecule is available. Used for mapq calculations only.")
var alignment_score_factor = flag.Float64("alignment_score_factor", 1.0, "Alignment score factor")
var genome = flag.String("genome", "", "Genome FASTA path")
var SIMULATED_DATA = flag.Bool("simulated", false, "Are data simulated?")
var output = flag.String("output", "", "Full path at which to output BAM file")
var read_group = flag.String("read_group", "default_read_group", "Read group")
var sample_id = flag.String("sample_id", "default_sample_id", "Sample name")
var threads = flag.Int("threads", 8, "Threads to use")
var max_bcs = flag.Int("max_bcs", -1, "Maximum nubmer of barcodes to process")
var DEBUG = flag.Bool("debug", false, "Debug mode")
var positionChunkSize = flag.Int("position_chunk_size", 40000000, "Bases across which to chunk within a chromosome for the purposes of bucketing by barcode, sorting, merging, so that we can do a fast samtools cat on the final BAMs.")
var debugTags = flag.Bool("debugBamTags", false, "Debug BAM tags")

type ChainedHit struct {
	contig    string
	pos       int64
	aend      int64
	read_id   int
	mate_id   int
	hit_id    int
	secondary bool
	read1     bool
	score     int
	chain     unsafe.Pointer //*C.mem_alnreg_t
	read      *[]byte
	fastq     *fastqreader.FastQRecord
	aln       *EasyAlignment
}

type Alignment struct {
	id                                int
	read1                             bool
	is_proper                         bool
	soft_clipped                      bool
	raw_barcode                       *[]byte
	barcode_qual                      *[]byte
	read_name                         *string
	read_seq                          *[]byte
	read_qual                         *[]byte
	sample_index                      *[]byte
	sample_index_qual                 *[]byte
	mapq                              int
	molecule_difference               float64
	contig                            string
	pos                               int64
	aend                              int64
	score                             int
	mismatches                        int
	matches                           int
	indels                            int
	read_id                           int
	bad_molecule                      bool
	correctly_placed                  bool
	mate_id                           int
	mate_alignment                    *Alignment
	reversed                          bool
	molecule_id                       int
	cigar                             []uint32
	read_group                        *string
	active                            bool    // the selected alignment for this read
	log_alignment_probability         float64 // does not include penalty for improperly paired
	updated_log_alignment_probability float64
	bwa_pick                          bool
	mapq_data                         *MapQData
	sum_move_probability_change       float64
	molecule_confidence               float64
	active_molecule                   bool
	readmap_s                         int
	readmap_e                         int
	secondary                         *Alignment
	duplicate                         bool
}

func (aln *Alignment) IsUnmapped() bool {
	if !aln.is_proper && aln.score-17 < 19 {
		return true
	}
	return false
}

type MapQData struct {
	copies                          int
	copies_in_active_molecules      int
	unique_molecules_active         int
	copies_outside_active_molecules int
	reads_in_molecule               int
	active_alignments_in_molecules  string
	second_best                     *Alignment
	second_best_proper_pair         bool
	second_best_molecule_reads      int
	second_best_molecule_confidence float64
}

/*
 Holds configuration parameters for the optimizer. These to be constant
 after set in main.
*/
type RFAConfig struct {
	improper_penalty float64
}

// types and functions to be able to sort a list of aligntments by position,
// must already be of the same contig otherwise it doesnt make sense
type ByPosition []*Alignment

func (a ByPosition) Len() int           { return len(a) }
func (a ByPosition) Swap(i, j int)      { a[i], a[j] = a[j], a[i] }
func (a ByPosition) Less(i, j int) bool { return a[i].pos < a[j].pos }

type CandidateMolecule struct {
	id                      int
	chrom                   string
	start                   int64
	stop                    int64
	alignments              *OrderedMap          //map[int]map[int]*Alignment read id to map of alignment id to alignment
	best_alignment_for_read *OrderedAlignmentMap //map[int]*Alignment read id to alignment
	active_alignments       *OrderedAlignmentMap //map[int]*Alignment read id to alignment
	log_probability         float64
	true_molecule           bool
	molecule_confidence     float64
	differences             float64
	soft_clipped            int
}

type Optimizer struct {
	candidate_molecules       []*CandidateMolecule
	alignments                [][]*Alignment
	total_alignment_score     float64
	currentMoleculeMoveSource int
	currentScore              float64
	log_molecule_probability  float64
	log_unpaired_probability  float64
	barcode                   string
}

/*
 * Holds a single unit of "work" to be processed by an RFA thread
 */
type WorkUnit struct {
	reads          []fastqreader.FastQRecord
	barcodenum     int
	unique_barcode bool
}

/*
 * Holds statistics and coordination points for RFA
 */
type RFAStats struct {
	total                 int64
	correct               int64
	correct_mapq10        int64
	total_mapq10          int64
	total_improper        int64
	total_improper_before int64
	mapq                  int64

	/* This lock synchronizes access to the "mapq.csv" file */
	file *bufio.Writer
	lock sync.Mutex
}

var barcode_reads [][]fastqreader.FastQRecord
var barcode_reads_lock sync.Mutex

func ReturnBuffer(reads []fastqreader.FastQRecord) {
	barcode_reads_lock.Lock()
	barcode_reads = append(barcode_reads, reads)
	barcode_reads_lock.Unlock()
}

type Data struct {
	alignments [][]*Alignment
	reads      []fastqreader.FastQRecord
	attach_bx  bool
}

//TODO remove TODO's when fixing thing TODO's relate to, not long after lol
func main() {

	// Use worker thread count request on cmdline, or
	// all CPUs if -threads wasn't specified
	numCPU := runtime.NumCPU()
	if *threads < 0 {
		*threads = numCPU
	}
	runtime.GOMAXPROCS(*threads + 2)

	flag.Parse()

	if _, err := os.Stat(*reads); os.IsNotExist(err) {
		panic(fmt.Sprintf("File does not exist %s", *reads))
	}
	if _, err := os.Stat(*genome); os.IsNotExist(err) {
		panic(fmt.Sprintf("Fasta file not found %s", *genome))
	}
	if syscall.Access(*output, 2) != nil { //is output writable
		panic(fmt.Sprintf("Output directory not writable by this process %s", *output))
	}
	print(fmt.Sprintf("Starting lariat. Version: %s\n", __VERSION__))

	fastq, err := fastqreader.OpenFastQ(*reads)

	if err != nil {
		panic(err)
	}
	print(fmt.Sprintf("Loading reference genome: %s\n", *genome))

	ref := GoBwaLoadReference(*genome)
	print("Reference loaded\n")
	settings := GoBwaAllocSettings()
	config := &RFAConfig{}

	config.improper_penalty = float64(*improper_pair_penalty) * (*alignment_score_factor)

	var w *bufio.Writer

	barcode_num := 0
	bams, err := CreateBAMs(ref, *output, *read_group, *sample_id, *positionChunkSize, *debugTags)
	if err != nil {
		panic(err)
	}
	work_to_do := make(chan *WorkUnit, 2)
	//finished := make (chan bool);

	stats := &RFAStats{}
	stats.file = w

	/* This RW lock gets Rlocked by each worker thread. Worker threads
	 * unlock it when they are finished and have copied their data to the
	 * BAM writer
	 */
	var worker_lock sync.RWMutex

	barcode_reads = [][]fastqreader.FastQRecord{}

	/* Start some workers */
	for i := 0; i < *threads; i++ {
		go WorkerThread(work_to_do, bams, ref, settings, config, stats, &worker_lock)
	}

	/* Iterate over source file, giving work to the workers */
	for {
		barcode_num++
		if len(barcode_reads) == 0 {
			barcode_reads_lock.Lock()
			barcode_reads = append(barcode_reads, make([]fastqreader.FastQRecord, 0, 50000))
			barcode_reads_lock.Unlock()
		}
		barcode_reads_lock.Lock()
		bc_reads := barcode_reads[len(barcode_reads)-1]
		//fmt.Println(len(barcode_reads))
		barcode_reads = barcode_reads[0 : len(barcode_reads)-1]
		barcode_reads_lock.Unlock()
		bc_reads, err, full_barcode := fastq.ReadBarcodeSet(&bc_reads)
		if err != nil {
			break
		}

		// Max BCs to process
		if barcode_num == *max_bcs {
			break
		}

		work_to_do <- &WorkUnit{bc_reads, barcode_num, full_barcode}
	}

	/* Tell each worker to exit */
	for i := 0; i < *threads; i++ {
		work_to_do <- nil
	}

	/* Wait for each worker to finish any final tasks and exit */
	worker_lock.Lock()
	worker_lock.Unlock()

	/* Close and flush the BAM file */
	bams.Close()
	fmt.Println("Lariat completed successfully")
}

/*
 * This is a single "worker" thread. It tries to grab work units until it gets
 * nil, then it shuts down.
 */
func WorkerThread(input chan *WorkUnit,
	bams *BAMWriters,
	ref *GoBwaReference,
	settings *GoBwaSettings,
	config *RFAConfig,
	stats *RFAStats,
	worker_lock *sync.RWMutex) {

	worker_lock.RLock()

	for work := <-input; work != nil; work = <-input {
		DoRFAForOneBarcode(work, bams, ref, settings, config, stats, work.reads)
	}
	worker_lock.RUnlock()
}

func StashAlignments(a *OrderedMap) *OrderedMap {

	res := NewOrderedMap()

	for id, read := range a.Iter() {

		a_for_read := NewOrderedMap()
		res.Set(id, a_for_read)

		alignments := FixGetForTypeOrderedMap(read)
		for id2, a_ptr := range alignments.Iter() {
			aln := FixGetForTypeAlignment(a_ptr)
			a_for_read.Set(id2, aln)
		}
	}

	return res
}

func DoRFAForOneBarcode(work *WorkUnit,
	bams *BAMWriters,
	ref *GoBwaReference,
	settings *GoBwaSettings,
	config *RFAConfig,
	stats *RFAStats,
	reads []fastqreader.FastQRecord) {

	stats.total = 0
	stats.mapq = 0
	barcode_num := work.barcodenum
	barcode_reads := work.reads
	arena := NewArena()
	worthRunningRFA := worthRunningRFA(barcode_reads, work.unique_barcode)
	barcode_chains, barcode := getChains(ref, settings, barcode_reads, arena, 25, 17)
	alignments, stashed_alignments := getAlignments(ref, settings, barcode_chains, *alignment_score_factor, 17, arena)
	//stashed_alignments := StashAlignments(alignments);

	positions := tagBestAlignments(alignments, -17)

	//print("considering working on barcode " + string(barcode_reads[0].Barcode10X) + "\n")
	if !worthRunningRFA {
		estimateMapQualities(-1, alignments, nil, config.improper_penalty, 0, stats)
		markDuplicates(alignments)
		CheckSplitReads(stashed_alignments)
		DumpToBams(&Data{alignments: alignments, reads: reads, attach_bx: work.unique_barcode}, bams)
		arena.Free()
		return
	}

	print("working on barcode " + string(barcode_reads[0].Barcode10X) + "\n")

	candidate_molecules := inferMolecules(positions)
	markBestAlignmentForReadInMolecule(candidate_molecules, config.improper_penalty)
	candidate_molecules = scrapMolecules(candidate_molecules)

	setMoleculeDifferences(candidate_molecules, false)

	optimizer_obj := &Optimizer{
		candidate_molecules:       candidate_molecules,
		alignments:                alignments,
		currentMoleculeMoveSource: 0,
		log_molecule_probability:  *molecule_penalty,
		log_unpaired_probability:  config.improper_penalty,
		barcode:                   barcode,
	}

	optimizer_obj.currentScore = globalScore(optimizer_obj)

	optimized := optimizer.Optimize(optimizer.Optimizable(*optimizer_obj), 1, 2, 4*len(candidate_molecules)).(Optimizer)

	estimateMapQualities(barcode_num, optimized.alignments, optimized.candidate_molecules, optimized.log_unpaired_probability, optimized.log_molecule_probability, stats)

	if *SIMULATED_DATA {
		//check
		for j := 0; j < len(optimized.candidate_molecules); j++ {
			molecule := optimized.candidate_molecules[j]
			for _, alignment := range molecule.active_alignments.Iter() {
				atomic.AddInt64(&(stats.total), 1)
				if alignment.mapq >= 10 {
					atomic.AddInt64(&(stats.total_mapq10), 1)
				}
				var pos float64
				if alignment.read1 {
					pos, _ = strconv.ParseFloat(strings.Split(*alignment.read_name, ":")[5], 64)
				} else {
					pos, _ = strconv.ParseFloat(strings.TrimSpace(strings.Split(*alignment.read_name, ":")[6]), 64)
				}

				if math.Abs(pos-float64(alignment.pos)) < 600 {
					atomic.AddInt64(&(stats.correct), 1)
					alignment.correctly_placed = true
					if alignment.mapq >= 10 {
						atomic.AddInt64(&(stats.correct_mapq10), 1)
					}
				}
			}
		}
	}
	markDuplicates(alignments)
	CheckSplitReads(stashed_alignments)
	DumpToBams(&Data{optimized.alignments, reads, true}, bams)
	arena.Free()
}

func DeAlignCrappyReads(reads [][]*Alignment) {
	for _, readArray := range reads {
		for _, aln := range readArray {
			if !aln.is_proper && aln.score-17 < 19 {
				aln.pos = -1
			}
		}
	}
}

func unbarcodeAlignments(alignments [][]*Alignment) {
	for _, alignmentList := range alignments {
		for _, alignment := range alignmentList {
			if alignment.active {
				barcode := []byte(strings.Split(string(*alignment.raw_barcode), "-")[0])
				alignment.raw_barcode = &barcode
			}
		}
	}
}

func setMoleculeDifferences(candidate_molecules []*CandidateMolecule, setBad bool) {
	for i := 0; i < len(candidate_molecules); i++ {
		differences := 0
		for _, alignment := range candidate_molecules[i].active_alignments.Iter() {
			differences += alignment.mismatches
		}
		if setBad {
			if float64(differences)/float64(candidate_molecules[i].active_alignments.Len()) > 1.5 {
				for _, alignment := range candidate_molecules[i].active_alignments.Iter() {
					alignment.bad_molecule = true
				}
			}
		}
		candidate_molecules[i].differences = float64(differences) / float64(candidate_molecules[i].active_alignments.Len())
		for _, alignment := range candidate_molecules[i].active_alignments.Iter() {
			alignment.molecule_difference = candidate_molecules[i].differences
		}
	}

}

func scoreAlignment(aln *Alignment, mate *Alignment) float64 {
	score := float64(aln.mismatches*-2.0 + aln.indels*-3.0)
	if aln.soft_clipped {
		score += -5.0
	}
	if mate != nil {
		score += float64(mate.mismatches*-2.0 + mate.indels*-3.0)
		if mate.soft_clipped {
			score += -5.0
		}
	}
	if mate == nil || !isPair(aln, mate) {
		score += float64(*improper_pair_penalty) / 2.5
	}
	if !aln.active_molecule {
		score += *log10_inactive_molecule
	}
	return score
}

// If two reads have the same value, then they are duplicates
type readDupTuple struct {
	read1      bool
	reversed   bool
	contig     string
	pos        int64
	mateContig string
	matePos    int64
}

// For each read, make a tuple of (bc_sequence, read.is_read1, read.is_reverse, read.tid, read.pos, read.mrnm, read.mpos)
// reads with an equal value of this tuple are defined as duplicates or one another.
// mark all but 1 read in each group as a duplicate
func markDuplicates(alignments [][]*Alignment) {

	dupSeen := make(map[readDupTuple]bool)

	//now go through every read_id and normalize all alternate alignment probabilities
	for _, alignmentArray := range alignments {
		for _, alignment := range alignmentArray {
			if alignment.active {

				mateAlignment := alignment.mate_alignment

				readTuple := readDupTuple{
					read1:      alignment.read1,
					reversed:   alignment.reversed,
					contig:     alignment.contig,
					pos:        alignment.pos,
					mateContig: mateAlignment.contig,
					matePos:    mateAlignment.pos}

				// If we have seen this tuple before, mark it as duplicate
				// Otherwise note tuple
				_, haveSeen := dupSeen[readTuple]
				if haveSeen {
					alignment.duplicate = true
				} else {
					dupSeen[readTuple] = true
				}
			}
		}
	}
}

// Basic strategy is two-fold:
// 1. is based on calculating a mapq vs being able to move that read to any of the other
//    alignments that it has. This includes the probabilities of a read being a singleton vs being in an active molecule.
//    Every alignment has some probability calculated from its alignment score/cigar string. Then alignments are
//    penalized for improper_pair and for being outside of active molecules. The mapq is then calculated by normalizing the sum of
//    alignment probabilities to 1 and then the mapq is -10*log10(1-probability_of_chosen_alignment)
//    If there is a problem with this method it is that it dings alignments for not being in an active read, but if a whole
//    molecule could have moved to that location with very little probability change making it an active molecule, the read
//    shouldn't have gotten that ding. I have tried not dinging reads in candidate molecules that are inactive but have
//    greater than some number of alignments but have found it better to just include method #2
// 2. is a mapq considering the probability change in the case of a candidate molecule making a sub-move (moving all of the reads
//    that can move from this molecule to another candidate molecule) to every candidate molecule that it can make a sub move to.
//    So for each active molecule we calculate the probability change of a sub-move to every candidate molecule that it can. For each
//    sub-move calculated, every alignment that was involved in that sub move gets that probability change added to a sum.
//    The log probability change of making a sub-move from a molecule to itself is 0 so probability 1.0. then we normalize the probability
//    of the final selected sub-move which is the molecule moving to itself p = 1.0 vs the sum of all of the probabilities of the moves
//    which is p_final = 1.0/sum_i(p(move_i)) and then calculate a mapq from that probability = -10*log10(1-p_final)
//Finally we take the min of the two approaches as the final mapq
//
func estimateMapQualities(barcode int,
	alignments [][]*Alignment,
	candidate_molecules []*CandidateMolecule,
	log_unpaired_probability, log_molecule_probability float64,
	stats *RFAStats) {
	improper := 0
	proper := 0
	read_copies_in_active_molecule := map[int]int{}     //TODO remove, book keeping
	read_copies_not_in_active_molecule := map[int]int{} //TODO remove, book keeping
	unique_molecules_active := map[int]map[int]bool{}
	// calculate whether reads are proper pairs

	for read_id, alignmentArray := range alignments {
		for _, alignment := range alignmentArray {
			mate_id := alignment.mate_id
			is_proper := false
			mapqData := MapQData{active_alignments_in_molecules: ""}
			alignment.mapq_data = &mapqData
			mateAlignments := alignments[mate_id]
			for _, mateAlignment := range mateAlignments {
				if alignment.active && mateAlignment.active {
					alignment.mate_alignment = mateAlignment
				}
				if isPair(alignment, mateAlignment) {
					if alignment.active {
						if mateAlignment.active {
							is_proper = true
							proper++
						} else {
							improper++
						}
					} else {
						is_proper = true
					}
					if alignment.mate_alignment == nil {
						alignment.mate_alignment = mateAlignment
					}
				}
			}
			alignment.updated_log_alignment_probability = alignment.log_alignment_probability
			//ding alignment probability for being improper pair
			if !is_proper {
				alignment.updated_log_alignment_probability += (log_unpaired_probability)
			}
			alignment.is_proper = is_proper
		}
		read_copies_in_active_molecule[read_id] = 0
		read_copies_not_in_active_molecule[read_id] = 0
		unique_molecules_active[read_id] = map[int]bool{}
	}
	// Now to update alignment probabilities for being singleton/outside active molecules
	// this part only happens if we ran RFA, bad barcodes etc get no more probability updates
	if candidate_molecules != nil {
		setMoleculeConfidences(candidate_molecules, log_unpaired_probability, log_molecule_probability)
		setMoleculeDifferences(candidate_molecules, false)

		// #1 update alignment probabilities based on whether they are in an active molecule or not
		// here we take 2 sources of power. 1 is taht singletons are rare
		for read_id, alignmentArray := range alignments {
			for _, alignment := range alignmentArray {
				is_molecule_active := false
				if alignment.molecule_id != -1 {
					molecule := candidate_molecules[alignment.molecule_id]
					is_molecule_active = molecule.active_alignments.Len()-molecule.soft_clipped > 4 && molecule.molecule_confidence > 0.1
					alignment.active_molecule = is_molecule_active
				}
				if is_molecule_active {
					read_copies_in_active_molecule[read_id]++ //TODO remove, book keeping
					unique_molecules_active[read_id][alignment.molecule_id] = true
				} else {
					read_copies_not_in_active_molecule[read_id]++ //TODO remove, book keeping
				}
				if alignment.molecule_id != -1 {
					alignment.mapq_data.reads_in_molecule = candidate_molecules[alignment.molecule_id].active_alignments.Len()
				}
			}
		}
	}
	if *SIMULATED_DATA {
		for _, alignmentArray := range alignments {
			for _, alignment := range alignmentArray {
				name := strings.Split(*alignment.read_name, ":")
				pos := int64(0)
				if alignment.read1 {
					pos, _ = strconv.ParseInt(name[5], 10, 64)
				} else {
					blah, _ := strconv.ParseFloat(strings.TrimSpace(strings.Split(*alignment.read_name, ":")[6]), 64)
					pos = int64(blah)
				}
				chrom := name[2]
				for i := 0; i < len(candidate_molecules); i++ {
					mol := candidate_molecules[i]
					if mol.chrom == chrom {
						if mol.stop >= pos && mol.start <= pos {
							mol.true_molecule = true
						}
					}
				}
			}
		}
	}
	//now go through every read_id and normalize all alternate alignment probabilities
	for read_id, alignmentArray := range alignments {
		second_best_proper_pair := false
		second_best_log_probability := -1000.0
		second_best_molecule_reads := -1
		var second_best_alignment *Alignment
		second_best_molecule_confidence := -1.0
		total_probability := float64(0)
		scores := []float64{}
		for _, alignment := range alignmentArray {
			mateArray := alignments[alignment.mate_id]
			best_score := -math.MaxFloat64
			for _, mateAlignment := range mateArray {
				score := scoreAlignment(alignment, mateAlignment)
				if !alignment.active {
					if score > second_best_log_probability {
						second_best_log_probability = alignment.updated_log_alignment_probability
						second_best_alignment = alignment
						alignment.mate_alignment = mateAlignment
						second_best_proper_pair = alignment.is_proper
						if alignment.molecule_id != -1 {
							alt_mol := candidate_molecules[alignment.molecule_id]
							second_best_molecule_confidence = alt_mol.molecule_confidence
							second_best_molecule_reads = alt_mol.active_alignments.Len()
						}
					}
				}

				if score > best_score {
					best_score = score
					// fmt.Println(alignment.contig, alignment.pos, alignment.score,"mismatches", alignment.mismatches,"indels",alignment.indels,"soft",alignment.soft_clipped, alignment.cigar,
					// 	 mate_alignment.contig, mate_alignment.pos, mate_alignment.score, "mismatches",mate_alignment.mismatches,"indels",mate_alignment.indels, "soft",mate_alignment.soft_clipped,mate_alignment.cigar,"proper_pair?" ,isPair(alignment, mate_alignment),"score",score)

				}

			}
			if len(mateArray) == 0 {
				best_score = scoreAlignment(alignment, nil)
			}
			scores = append(scores, best_score)
		}

		for _, alignment := range alignmentArray {
			if alignment.active {
				alignment.mapq_data.second_best = second_best_alignment
				alignment.mapq_data.second_best_proper_pair = second_best_proper_pair
				alignment.mapq_data.second_best_molecule_confidence = second_best_molecule_confidence
				alignment.mapq_data.second_best_molecule_reads = second_best_molecule_reads
			}
		}

		sort.Float64s(scores)
		for i := len(scores) - 1; i >= 0 && len(scores)-i <= 15; i-- {
			total_probability += math.Pow(10, scores[i])
		}
		debug_strings := map[int]map[int]string{}
		for _, alignment := range alignmentArray {

			if alignment.active {
				alignment.mapq_data.second_best_proper_pair = second_best_proper_pair
				alignment.mapq_data.second_best_molecule_reads = second_best_molecule_reads
				alignment.mapq_data.second_best_molecule_confidence = second_best_molecule_confidence
				alignment.mapq_data.second_best = second_best_alignment
				alignment.mapq_data.copies = len(alignmentArray)
				alignment.mapq_data.copies_in_active_molecules = read_copies_in_active_molecule[alignment.read_id]
				alignment.mapq_data.copies_outside_active_molecules = read_copies_not_in_active_molecule[read_id]
				alignment.mapq_data.unique_molecules_active = len(unique_molecules_active[read_id])
				if *DEBUG {
					alt_alignments := alignments[alignment.read_id]
					for _, alignment_alt := range alt_alignments {
						if alignment_alt.molecule_id != -1 {
							chrom := alignment_alt.contig
							start := candidate_molecules[alignment_alt.molecule_id].start
							end := candidate_molecules[alignment_alt.molecule_id].stop
							sinksource := 0
							sourcesink := 0
							molstring := ""
							_, has := debug_strings[alignment.molecule_id]
							if has {
								_, has = debug_strings[alignment.molecule_id][alignment_alt.molecule_id]
							}
							if !has {
								for _, read1 := range candidate_molecules[alignment.molecule_id].active_alignments.Iter() {
									read := candidate_molecules[alignment_alt.molecule_id].best_alignment_for_read.Get(read1.id)
									if candidate_molecules[alignment.molecule_id].start == 30943998 && read != nil {
										fmt.Println(string(*read1.read_name), "source", read1.score, read1.contig, read1.pos, "sink", read.score, read.contig, read.pos)
									}
									if read != nil {
										sourcesink++
									}
								}
								for _, rid := range candidate_molecules[alignment_alt.molecule_id].active_alignments.IterKeys() {
									has := FixGetForTypeAlignment(candidate_molecules[alignment.molecule_id].best_alignment_for_read.Get(rid)) != nil
									// read := FixGetForTypeAlignment(candidate_molecules[alignment.molecule_id].best_alignment_for_read.Get(rid))
									// read_right := FixGetForTypeAlignment(candidate_molecules[alignment_alt.molecule_id].best_alignment_for_read.Get(rid))

									if has {
										sinksource++
									}
								}

								ST := strconv.FormatInt(int64(sourcesink), 10)
								TS := strconv.FormatInt(int64(sinksource), 10)
								sourcesinkchange, _ := fastScore(candidate_molecules[alignment.molecule_id], candidate_molecules[alignment_alt.molecule_id], log_unpaired_probability, log_molecule_probability, false)

								sinksourcechange, _ := fastScore(candidate_molecules[alignment_alt.molecule_id], candidate_molecules[alignment.molecule_id], log_unpaired_probability, log_molecule_probability, false)
								active := strconv.FormatInt(int64(candidate_molecules[alignment_alt.molecule_id].active_alignments.Len()), 10)
								spots := strconv.FormatInt(int64(candidate_molecules[alignment_alt.molecule_id].best_alignment_for_read.Len()), 10)
								STC := strconv.FormatInt(int64(sourcesinkchange), 10)
								TSC := strconv.FormatInt(int64(sinksourcechange), 10)
								molstring = chrom + ":" + strconv.FormatInt(start, 10) + "-" + strconv.FormatInt(end, 10) + ":alignments:" + active + ":spots:" + spots + ":mv_S->T:" + ST + ":" + STC + ":mv_T->S:" + TS + ":" + TSC + ","
								_, has_key := debug_strings[alignment.molecule_id]
								if !has_key {
									debug_strings[alignment.molecule_id] = map[int]string{}
								}
								debug_strings[alignment.molecule_id][alignment_alt.molecule_id] = molstring
							} else {
								molstring = debug_strings[alignment.molecule_id][alignment_alt.molecule_id]
							}
							alignment.mapq_data.active_alignments_in_molecules += molstring

						}
					}
				}
			}
			score := scoreAlignment(alignment, alignment.mate_alignment)
			mapq := -10.0 * math.Log10(1.0-math.Pow(10, score)/total_probability)

			mapq = math.Min(float64(60), mapq)

			alignment.mapq = int(mapq)
		}
	}
}

func setMoleculeConfidences(molecules []*CandidateMolecule, log_unpaired_probability, log_molecule_probability float64) {
	for i := 0; i < len(molecules); i++ {
		molecules[i].molecule_confidence = moleculeConfidence(molecules[i], molecules[i].active_alignments.Len())
		for _, alignment := range molecules[i].active_alignments.Iter() {

			if alignment.soft_clipped {
				molecules[i].soft_clipped++
			}
			alignment.molecule_confidence = molecules[i].molecule_confidence
		}
	}
}

func scrapMolecules(candidate_molecules []*CandidateMolecule) []*CandidateMolecule {
	toReturn := []*CandidateMolecule{}
	count := 0
	for i := 0; i < len(candidate_molecules); i++ {
		if candidate_molecules[i].active_alignments.Len() > 0 {
			toReturn = append(toReturn, candidate_molecules[i])
			for _, read_id := range candidate_molecules[i].alignments.IterKeys() {
				alignment_map := FixGetForTypeOrderedMap(candidate_molecules[i].alignments.Get(read_id))
				for _, alignment_id := range alignment_map.IterKeys() {
					alignment := FixGetForTypeAlignment(alignment_map.Get(alignment_id))
					alignment.molecule_id = count
				}
			}
			count++
		} else {
			for _, read_id := range candidate_molecules[i].alignments.IterKeys() {
				alignment_map := FixGetForTypeOrderedMap(candidate_molecules[i].alignments.Get(read_id))
				for _, alignment_id := range alignment_map.IterKeys() {
					alignment := FixGetForTypeAlignment(alignment_map.Get(alignment_id))
					alignment.molecule_id = -1
				}
			}
		}
	}
	return toReturn
}

func worthRunningRFA(barcode_reads []fastqreader.FastQRecord, uniqueBarcode bool) bool {
	if len(barcode_reads) == 0 || !uniqueBarcode {
		return false
	}
	bcParts := strings.Split(string(barcode_reads[0].Barcode10X), "-")
	if len(bcParts) < 2 {
		return false
	}
	if len(barcode_reads) < 50 {
		return false
	}
	return true
}

func isPair(read1, read2 *Alignment) bool {
	if read1.reversed == read2.reversed || read1.contig != read2.contig {
		return false
	}
	var forward, reverse *Alignment
	if read1.reversed {
		forward = read2
		reverse = read1
	} else {
		forward = read1
		reverse = read2
	}
	dist := reverse.pos - forward.pos
	//TODO dont delete, trimming code to be turned on at later date
	// if dist < 0 && dist >= -35 {
	// 	if reverse.cigar[0] == uint32(3) && len(reverse.cigar) > 2 && reverse.cigar[2] == 0 && int64(reverse.cigar[3]) > -dist {
	// 		reverse.cigar[1] += uint32(-dist)
	// 		reverse.cigar[3] -= uint32(-dist)
	// 		reverse.pos += -dist
	// 	}
	// }
	// 	cigar_end := len(forward.cigar)
	// 	overhang := reverse.aend - forward.aend
	// 	if overhang < 0 {
	// 		if forward.cigar[cigar_end-2] == uint32(3) && cigar_end > 2 && int64(forward.cigar[cigar_end-3]) > -overhang {
	// 			forward.cigar[cigar_end-1] += uint32(-overhang)
	// 			forward.cigar[cigar_end-3] -= uint32(-overhang)
	// 			forward.aend -= -overhang
	// 		}
	// 	}
	return dist >= int64(-35) && dist < int64(750)
}

func (o Optimizer) GenerateMove(accept_move func(p_curr float64, p_next float64) bool) optimizer.Optimizable {
	sourceMolecule := o.candidate_molecules[o.currentMoleculeMoveSource]

	if sourceMolecule.active_alignments.Len() == 0 {
		o.currentMoleculeMoveSource = (o.currentMoleculeMoveSource + 1) % len(o.candidate_molecules)
		return o
	}
	var sinkMolecule *CandidateMolecule
	best_move := Move{score_change: -math.MaxFloat64}

	for i := 0; i < len(o.candidate_molecules); i++ {
		if i == o.currentMoleculeMoveSource {
			continue
		}
		sinkMolecule = o.candidate_molecules[i]

		score, move := o.fastScore(sourceMolecule, sinkMolecule, false)

		if (score > best_move.score_change ||
			(score == best_move.score_change && move.sink.active_alignments.Len() > best_move.sink.active_alignments.Len())) && move.num_moved > 0 {
			best_move = move
		}
	}

	best_score := best_move.score_change
	if best_move.alignment_change > 30 {
		_, best_move = o.fastScore(sourceMolecule, best_move.sink, true)
	}

	if best_score > 0 || (best_score == 0 && best_move.sink.active_alignments.Len() > sourceMolecule.active_alignments.Len()) {
		acceptMove(best_move)
	}
	//newScoreGlobal := globalScore(&o)
	//newScoreFast := o.currentScore + best_move.score_change

	//scoreCmp := fmt.Sprintf("Global Score: %f  FastScore:%f", newScoreGlobal, newScoreFast)
	//print(scoreCmp)

	o.currentMoleculeMoveSource = (o.currentMoleculeMoveSource + 1) % len(o.candidate_molecules)
	return o
}

type Move struct {
	score_change     float64
	alignment_change float64
	source           *CandidateMolecule
	sink             *CandidateMolecule
	toDelete         []int
	toSet            []*Alignment
	num_moved        int
}

func globalScore(o *Optimizer) float64 {
	singleton_penalty := 0.01
	score := 0.0

	for _, molecule := range o.candidate_molecules {
		confidence := moleculeConfidence(molecule, molecule.active_alignments.Len())

		mol_score := o.log_molecule_probability * confidence
		score += mol_score

		for _, _al := range molecule.active_alignments.Iter() {
			al := FixGetForTypeAlignment(_al)
			alP := math.Exp(al.log_alignment_probability)
			score += math.Log(confidence*alP + (1-confidence)*alP*singleton_penalty)

			mate := FixGetForTypeAlignment(molecule.active_alignments.Get(al.mate_id))

			if mate == nil || !isPair(al, mate) {
				score += o.log_unpaired_probability
			}
		}
	}

	return score
}

func fastScore(sourceMolecule, sinkMolecule *CandidateMolecule, log_unpaired_probability, log_molecule_probability float64, strict bool) (float64, Move) {
	change := float64(0)
	alignment_change := float64(0)
	num := 0
	toDelete := []int{}
	toSet := []*Alignment{}
	soft_clipped := 0
	for _, sourceAlignment := range sourceMolecule.active_alignments.Iter() {
		if sourceAlignment.soft_clipped {
			soft_clipped += 1
		}
		read_id := sourceAlignment.read_id
		sinkAlignment := sinkMolecule.best_alignment_for_read.Get(read_id)
		if sinkAlignment != nil {

			mate_id := sourceAlignment.mate_id

			sourceMate := sourceMolecule.active_alignments.Get(mate_id)
			source_has_mate := sourceMate != nil
			source_has_mate_pair := source_has_mate && isPair(sourceAlignment, sourceMate)
			mate := sinkMolecule.best_alignment_for_read.Get(mate_id)
			sink_has_mate_pair := mate != nil && isPair(sinkAlignment, mate) && source_has_mate
			if strict && sourceAlignment.log_alignment_probability > sinkAlignment.log_alignment_probability {
				continue
			}

			if !source_has_mate_pair || (source_has_mate && sink_has_mate_pair) {
				if strict && source_has_mate_pair && mate.log_alignment_probability < sourceMate.log_alignment_probability {
					continue
				}
				toDelete = append(toDelete, read_id)
				toSet = append(toSet, sinkAlignment)
			}
			alignment_change += sinkAlignment.log_alignment_probability - sourceAlignment.log_alignment_probability
			if source_has_mate_pair && !sink_has_mate_pair && sourceMolecule.id != sinkMolecule.id {

				alignment_change += log_unpaired_probability / 2.0
			} else if !source_has_mate_pair && sink_has_mate_pair && sourceMolecule.id != sinkMolecule.id {
				alignment_change -= log_unpaired_probability / 2.0
			}
			num++
		}
	}

	source_active_before := isActiveMolecule(sourceMolecule, 0)
	source_active_after := isActiveMolecule(sourceMolecule, -num)
	if !source_active_after && source_active_before && sourceMolecule.id != sinkMolecule.id {
		change -= float64(sourceMolecule.best_alignment_for_read.Len()) * -0.5
	}
	sink_active_before := isActiveMolecule(sinkMolecule, 0)
	sink_active_after := isActiveMolecule(sinkMolecule, num)
	if sink_active_after && !sink_active_before && sourceMolecule.id != sinkMolecule.id {
		change += float64(sinkMolecule.best_alignment_for_read.Len()) * -0.5
	}
	if sourceMolecule.active_alignments.Len()-num == 0 && num > 0 && sourceMolecule.id != sinkMolecule.id {
		change -= -3.0
	}
	if sinkMolecule.active_alignments.Len() == 0 && num > 0 && sourceMolecule.id != sinkMolecule.id {
		change += -3.0
	}
	change += alignment_change

	return change, Move{source: sourceMolecule, sink: sinkMolecule, toDelete: toDelete, toSet: toSet, num_moved: num, score_change: change, alignment_change: alignment_change}
}

func isActiveMolecule(mol *CandidateMolecule, read_change int) bool {
	active := float64(mol.active_alignments.Len() + read_change)
	potential := float64(mol.best_alignment_for_read.Len())
	if active <= 4 {
		return false
	}
	if active/potential < 0.1 {
		return false
	}
	return true
}

func (o Optimizer) fastScore(sourceMolecule, sinkMolecule *CandidateMolecule, strict bool) (float64, Move) {
	change, move := fastScore(sourceMolecule, sinkMolecule, o.log_unpaired_probability, o.log_molecule_probability, strict)
	return change, move
}

func moleculeConfidence(mol *CandidateMolecule, num_active int) float64 {
	density := float64(num_active) / float64(mol.best_alignment_for_read.Len())
	return density
}

func acceptMove(move Move) {
	toDelete := move.toDelete
	toSet := move.toSet
	for i := 0; i < len(toDelete); i++ {
		read_id := toDelete[i]
		sinkAlignment := toSet[i]
		sourceAlignment := move.source.active_alignments.Get(read_id)
		move.source.active_alignments.Delete(read_id)
		move.sink.active_alignments.Set(read_id, sinkAlignment)
		sourceAlignment.active = false
		sinkAlignment.active = true
	}
}

func inferMolecules(positions [][]*Alignment) []*CandidateMolecule {
	toReturn := []*CandidateMolecule{}
	molecule_num := 0
	var currentMolecule *CandidateMolecule
	for _, position_list := range positions {
		for i := 0; i < len(position_list); i++ {
			if i == 0 || (i > 0 && position_list[i].pos-position_list[i-1].pos > 50000) {
				if i > 0 {
					currentMolecule.stop = position_list[i-1].pos
				}
				currentMolecule = &CandidateMolecule{
					chrom:               position_list[i].contig,
					start:               position_list[i].pos,
					id:                  molecule_num,
					alignments:          NewOrderedMap(),
					molecule_confidence: 1.0,
				}
				aln_map := NewOrderedMap()
				aln_map.Set(position_list[i].id, position_list[i])
				currentMolecule.alignments.Set(position_list[i].read_id, aln_map)

				toReturn = append(toReturn, currentMolecule)
				molecule_num++
			}
			alignment_map := FixGetForTypeOrderedMap(currentMolecule.alignments.Get(position_list[i].read_id))

			if alignment_map != nil {
				alignment_map.Set(position_list[i].id, position_list[i])
			} else {
				aln_map := NewOrderedMap()
				aln_map.Set(position_list[i].id, position_list[i])
				currentMolecule.alignments.Set(position_list[i].read_id, aln_map)
			}
		}
		if len(position_list) > 0 {
			currentMolecule.stop = position_list[len(position_list)-1].pos
		}
	}
	return toReturn
}

func markBestAlignmentForReadInMolecule(molecules []*CandidateMolecule, improper_pair_penalty float64) {
	active_alignment_num := 0
	for i := 0; i < len(molecules); i++ {
		molecule := molecules[i]
		active_alignments := NewOrderedAlignmentMap()
		best_alignment_for_read := NewOrderedAlignmentMap()
		for _, read_id := range molecule.alignments.IterKeys() {
			alignments := FixGetForTypeOrderedMap(molecule.alignments.Get(read_id))
			best_score := -math.MaxFloat64
			var best_alignment *Alignment
			for _, alignment_id := range alignments.IterKeys() {
				alignment := FixGetForTypeAlignment(alignments.Get(alignment_id))
				mate_alignments := FixGetForTypeOrderedMap(molecule.alignments.Get(alignment.mate_id))
				if mate_alignments != nil && mate_alignments.Len() > 0 {
					score := float64(0)
					// loop over mates to get best pair
					for _, mate_alignment_id := range mate_alignments.IterKeys() {
						mate_alignment := FixGetForTypeAlignment(mate_alignments.Get(mate_alignment_id))
						if isPair(alignment, mate_alignment) {
							score = alignment.log_alignment_probability + mate_alignment.log_alignment_probability
						} else {
							score = alignment.log_alignment_probability + mate_alignment.log_alignment_probability + improper_pair_penalty
						}
						if score > best_score {
							best_score = score
							best_alignment = alignment
						}
					}
				} else {
					if alignment.log_alignment_probability > best_score {
						best_score = alignment.log_alignment_probability
						best_alignment = alignment
					}
				}
				if alignment.active {
					active_alignments.Set(read_id, alignment)
					active_alignment_num++
				}
			}
			if best_alignment.active {
				active_alignments.Set(read_id, best_alignment)
			}
			best_alignment_for_read.Set(read_id, best_alignment)
		}
		molecule.active_alignments = active_alignments
		molecule.best_alignment_for_read = best_alignment_for_read
	}
}

// alignments sent back sorted by position
func tagBestAlignments(alignments [][]*Alignment, improper_pair_penalty float64) [][]*Alignment {
	positions := [][]*Alignment{}
	contigs := map[string]int{}
	read_ids_touched := make([]bool, len(alignments))
	improper_pair := 0
	proper_pair := 0
	reads := 0

	for read_id := range alignments {
		alignmentArray := alignments[read_id]
		touched := read_ids_touched[read_id]
		bestScore := -math.MaxFloat64
		first := true
		var bestAlignment *Alignment
		var bestMate *Alignment
		seed := int64(1)
		if len(alignmentArray) > 0 {
			md5sum := md5.Sum([]byte(*alignmentArray[0].read_name))
			seed = int64(binary.LittleEndian.Uint64(md5sum[0:8]))
		}
		rand.Seed(seed)
		reads++
		for _, alignment := range alignmentArray {
			if first {
				first = false
			}
			if read_id != alignment.read_id {
				panic(fmt.Sprintf("ids not idey%d%d", read_id, alignment.read_id))
			}
			mateAlignments := alignments[alignment.mate_id]
			totalScore := float64(0)
			for _, mateAlignment := range mateAlignments {
				//must consider the option of them separate vs them being together
				if isPair(alignment, mateAlignment) {
					random := (rand.Float64() / 2.0)
					totalScore = float64(alignment.score+mateAlignment.score) + random
				} else {
					totalScore = float64(alignment.score+mateAlignment.score) + improper_pair_penalty + (rand.Float64() / 2.0)
				}

				// fmt.Println("improper pair penalty",improper_pair_penalty, alignment.contig, alignment.pos, "mismatches", alignment.mismatches, "indels", alignment.indels, alignment.cigar, alignment.score, alignment.log_alignment_probability,
				// 	mate_alignment.contig, mate_alignment.pos, "paired?",isPair(alignment, mate_alignment), "mismatches",mate_alignment.mismatches, "indels",mate_alignment.indels,mate_alignment.cigar,mate_alignment.score, mate_alignment.log_alignment_probability,"\t", total_score)
				// fmt.Println(total_score,"pos", alignment.pos, mate_alignment.pos, "score", alignment.score, mate_alignment.score, "logp",alignment.log_alignment_probability, mate_alignment.log_alignment_probability, isPair(alignment, mate_alignment))
				if alignment.mate_id != mateAlignment.read_id {
					panic(fmt.Sprintf("mates not matey%d%d", alignment.mate_id, mateAlignment.read_id))
				}
				if totalScore > bestScore {
					bestScore = totalScore
					bestAlignment = alignment
					bestMate = mateAlignment
				}
			}
			if len(mateAlignments) == 0 {
				score := float64(alignment.score) + rand.Float64()/2.0
				if score > bestScore {
					bestScore = score
					bestAlignment = alignment
				}
			}
			index, has_chrom := contigs[alignment.contig]

			if has_chrom {
				positions[index] = append(positions[index], alignment)
			} else {
				contigs[alignment.contig] = len(positions)
				positions = append(positions, []*Alignment{alignment})
			}
		}
		if !touched {
			bestAlignment.active = true
			bestAlignment.bwa_pick = true

			if bestMate != nil {
				if !isPair(bestAlignment, bestMate) {
					improper_pair += 2
				} else {
					proper_pair += 2
					bestAlignment.is_proper = true
					bestMate.is_proper = true
				}
				bestMate.active = true
				bestMate.bwa_pick = true
				read_ids_touched[bestMate.read_id] = true
			} else {
				improper_pair += 1
			}
		}
	}
	for _, position := range positions {
		sort.Sort(ByPosition(position))
	}
	return positions
}

//returns a map from read id to a map of
func getAlignments(ref *GoBwaReference, settings *GoBwaSettings, barcode_chains [][]ChainedHit, score_factor float64, delta int, arena *Arena) ([][]*Alignment, [][]*Alignment) {

	toReturn := make([][]*Alignment, len(barcode_chains))
	full := make([][]*Alignment, len(barcode_chains))
	for i := 0; i < len(barcode_chains); i++ {
		bestScore := 0
		for _, chain := range barcode_chains[i] {
			if chain.score > bestScore {
				bestScore = chain.score
			}
		}

		for j := 0; j < len(barcode_chains[i]); j++ {
			chain := barcode_chains[i][j]
			var alignment SingleReadAlignment
			if chain.chain != nil {
				alignment = GoBwaSmithWaterman(ref, settings, string(*(barcode_chains[i][j].read)), chain.chain, arena)
			} else {
				alignment = SingleReadAlignment{}
			}

			matches := 0
			indels := 0
			indel_length := 0
			soft_clipping := 0
			soft_clipping_num := 0
			soft_clipping_length := 0
			for j := 0; j < len(alignment.Cigar); j += 2 {
				if alignment.Cigar[j] == 0 {
					matches += int(alignment.Cigar[j+1])
				} else if alignment.Cigar[j] == 1 || alignment.Cigar[j] == 2 {
					indels += 1
					indel_length += int(alignment.Cigar[j+1])
				} else if alignment.Cigar[j] == 3 {
					soft_clipping += 1
					soft_clipping_num += 1
					soft_clipping_length += int(alignment.Cigar[j+1])
				}
			}
			if soft_clipping > 0 {
				soft_clipping = 5
			}
			mismatches := alignment.EditDistance - indel_length
			matches -= mismatches
			if mismatches < 0 {
				mismatches = 0
			}
			log_probability_score := -2.0*float64(mismatches) + -3.0*float64(indels) + float64(-soft_clipping) - float64(soft_clipping_length)*0.5
			if soft_clipping_num > 1 {
				log_probability_score -= 5
			}

			var quals *[]byte
			if chain.read1 {
				quals = &chain.fastq.ReadQual1
			} else {
				quals = &chain.fastq.ReadQual2
			}
			pos := chain.pos
			aend := chain.aend
			if pos != -1 && alignment.Reversed {
				pos = chain.aend + 1
				aend = chain.pos + 1
			}
			full_alignment := Alignment{
				id:                          chain.hit_id,
				aend:                        aend,
				read_name:                   &chain.fastq.ReadInfo,
				read_seq:                    chain.read,
				read_qual:                   quals,
				matches:                     matches,
				mismatches:                  mismatches,
				indels:                      indels,
				soft_clipped:                soft_clipping > 0,
				read1:                       chain.read1,
				raw_barcode:                 &chain.fastq.Barcode10X,
				barcode_qual:                &chain.fastq.Barcode10XQual,
				contig:                      alignment.Chrom,
				pos:                         pos,
				molecule_id:                 -1,
				score:                       chain.score,
				cigar:                       alignment.Cigar,
				read_id:                     chain.read_id,
				mate_id:                     chain.mate_id,
				log_alignment_probability:   log_probability_score,
				reversed:                    alignment.Reversed,
				sample_index:                &chain.fastq.Barcode,
				sample_index_qual:           &chain.fastq.BarcodeQual,
				read_group:                  read_group,
				sum_move_probability_change: 1.0,
				molecule_confidence:         0.00075 * 0.025,
				duplicate:                   false,
			}

			if chain.aln != nil {
				full_alignment.readmap_s = chain.aln.ReadS
				full_alignment.readmap_e = chain.aln.ReadE
			}
			full[barcode_chains[i][j].read_id] = append(full[barcode_chains[i][j].read_id], &full_alignment)
			if full_alignment.score >= bestScore-delta {
				toReturn[barcode_chains[i][j].read_id] = append(toReturn[barcode_chains[i][j].read_id], &full_alignment)
			}
		}
	}
	return toReturn, full
}

func getChains(ref *GoBwaReference, settings *GoBwaSettings, reads_for_barcode []fastqreader.FastQRecord, arena *Arena, score_delta, improper_pair_penalty int) ([][]ChainedHit, string) {
	toReturn := [][]ChainedHit{}
	hit_num := 0
	var barcode string
	for i := 0; i < len(reads_for_barcode); i++ {
		if len(reads_for_barcode[i].Read1) == 0 || len(reads_for_barcode[i].Read2) == 0 {
			continue
		}
		read1_chains, read2_chains := GoBwaMemMateSW(ref, settings, &reads_for_barcode[i].Read1, &reads_for_barcode[i].Read2, arena, score_delta)
		barcode = string(reads_for_barcode[i].Barcode10X)
		read1_num := 0
		toReturn = append(toReturn, []ChainedHit{})
		for j := 0; j < len(read1_chains); j++ {
			read1_chain_n := ChainedHit{
				contig:    read1_chains[j].Contig,
				pos:       read1_chains[j].Offset,
				aend:      read1_chains[j].Alignment_end,
				read_id:   i * 2,
				mate_id:   i*2 + 1,
				hit_id:    hit_num,
				read1:     true,
				secondary: read1_chains[j].Secondary,
				score:     read1_chains[j].Score,
				chain:     unsafe.Pointer(read1_chains[j].ChainedHit),
				fastq:     &reads_for_barcode[i],
				read:      &reads_for_barcode[i].Read1,
				aln:       &read1_chains[j],
			}
			read1_num++
			toReturn[len(toReturn)-1] = append(toReturn[len(toReturn)-1], read1_chain_n)
			hit_num++
		}
		if read1_num == 0 {
			toReturn[len(toReturn)-1] = append(toReturn[len(toReturn)-1], ChainedHit{
				read_id: i * 2,
				mate_id: i*2 + 1,
				pos:     -1,
				read1:   true,
				chain:   nil,
				fastq:   &reads_for_barcode[i],
				read:    &reads_for_barcode[i].Read1,
			})
			hit_num++
		}
		toReturn = append(toReturn, []ChainedHit{})
		read2_num := 0
		for j := 0; j < len(read2_chains); j++ {
			read2_chain_n := ChainedHit{
				contig:    read2_chains[j].Contig,
				pos:       read2_chains[j].Offset,
				aend:      read2_chains[j].Alignment_end,
				read_id:   i*2 + 1,
				mate_id:   i * 2,
				hit_id:    hit_num,
				read1:     false,
				score:     read2_chains[j].Score,
				chain:     unsafe.Pointer(read2_chains[j].ChainedHit),
				secondary: read2_chains[j].Secondary,
				fastq:     &reads_for_barcode[i],
				read:      &reads_for_barcode[i].Read2,
				aln:       &read2_chains[j],
			}
			read2_num++
			toReturn[len(toReturn)-1] = append(toReturn[len(toReturn)-1], read2_chain_n)
			hit_num++
		}
		if read2_num == 0 {
			toReturn[len(toReturn)-1] = append(toReturn[len(toReturn)-1], ChainedHit{
				read_id: i*2 + 1,
				mate_id: i * 2,
				pos:     -1,
				hit_id:  hit_num,
				read1:   false,
				chain:   nil,
				fastq:   &reads_for_barcode[i],
				read:    &reads_for_barcode[i].Read2,
			})
			hit_num++
		}
	}
	return toReturn, barcode
}
