// Copyright (c) 2015 10X Genomics, Inc. All rights reserved.

package inference

import (
	"sort"
	// "log"
)

type SplitScoring struct {
	alignment *Alignment
	score     float64
}

type SortSplitScoring []SplitScoring

func (a SortSplitScoring) Len() int           { return len(a) }
func (a SortSplitScoring) Swap(i, j int)      { a[i], a[j] = a[j], a[i] }
func (a SortSplitScoring) Less(i, j int) bool { return a[i].score > a[j].score }

func abs(x int) int {
	if x < 0 {
		return -x
	} else {
		return x
	}
}

func GetSplitAlignment(primary *Alignment, alignments []*Alignment, centromeres map[string]Region) (*Alignment, float64) {

	//var normalizer float64;
	if primary.pos == -1 {
		return nil, 0.0
	}

	var Ps, Pe int

	/* Compute Ps, Pe the range of the read that was mapped */
	Ps = primary.readmap_s
	Pe = primary.readmap_e
	if Ps > Pe {
		Ps, Pe = Pe, Ps
	}

	/* Need at least 28 clipped bases to try to split */
	if (Pe - Ps) > len(*primary.read_seq)-15 {
		return nil, 0.0
	}
	/* Array to store split candidates */
	candidates := []SplitScoring{}

	/* Iterate over all of the other mappings of this read */
	for _, secondary_candidate := range alignments {

		/* Don't try to split a read with itsself*/
		if secondary_candidate.active {
			continue
		}

		/* Ignore unmapped stuff */
		if secondary_candidate.pos == -1 {
			continue
		}
		//Secondarystart, Secondaryend
		var Ss, Se int

		var overlap int

		/* Compute the part of the secondary read that is mapped */
		Ss = secondary_candidate.readmap_s
		Se = secondary_candidate.readmap_e
		if Ss > Se {
			Ss, Se = Se, Ss
		}

		if (Ps < Ss && Pe > Se) || (Ss < Ps && Se > Pe) {
			/* The secondary is fully contained in the primary (or visa versa)
			 * don't split.
			 */
			continue
		} else if Ps < Ss {
			overlap = int(Pe - Ss)
		} else {
			if Ps < Ss {
				panic("YEOWCH")
			}
			overlap = Se - Ps
		}

		/* If the overlap is small enough and if the secondary score is large enough
		 * consider this candidate */
		if overlap < (Se-Ss)/2 {
			//if (5 + secondary_candidate.score) > (primary.score-1)/2 {
			secondary_candidate.is_proper = isPair(secondary_candidate, primary.mate_alignment)
			if secondary_candidate.score >= 36 || secondary_candidate.is_proper {
				candidates = append(candidates, SplitScoring{secondary_candidate, float64(secondary_candidate.score)})
			}
		}
	}
	if len(candidates) == 0 {
		return nil, 0.0
	}

	/* Pick the candidate with the highest score */
	sort.Sort(SortSplitScoring(candidates))
	c := candidates[0].alignment

	var mapq float64

	/* Estimate the mapq score by comparing the best and second-best candidates */
    second_best := scoreAlignment(primary,nil,0.0) + psuedoCountAlignmentScore(candidates[0].alignment,0.0)
	if len(candidates) > 1 {
		mapq = float64(candidates[0].score - candidates[1].score)
        second_best = scoreAlignment(primary, candidates[1].alignment, 0.0)
	} else {
		mapq = float64(candidates[0].score)
	}

    centromereRegion, ok := centromeres[c.contig]
    start := -1
    end := -1
    if ok {
        start = centromereRegion.start
        end = centromereRegion.end
    }
    if c.pos > int64(start) && c.pos <= int64(end) {
        mapq = 0.0
    }

	if mapq > 60 {
		mapq = 60
	}

	c.mapq = int(mapq)

	return c, second_best
}

/*
 * Iterate over all reads and compute secondary "split" reads for some of them
 */
func CheckSplitReads(reads [][]*Alignment, centromeres map[string]Region) {
	for _, readArray := range reads {
		var active *Alignment
		for _, a := range readArray {
			if a.active {
				active = a
				break
			}
		}
		split, second_best := GetSplitAlignment(active, readArray, centromeres)
		active.secondary = split
        if split != nil {
            split.mapq_data = &MapQData{second_best_score: second_best, score: scoreAlignment(split, active.mate_alignment, 0.0) }
            split.primary = active
        }
	}
}
