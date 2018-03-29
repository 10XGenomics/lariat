// Copyright (c) 2016 10X Genomics, Inc. All rights reserved.

package test

import (
	"fastqreader"
	"gobwa"
	"inference"
	"testing"
)

func TestLariatZeroLengthRead(t *testing.T) {
	ref := gobwa.GoBwaLoadReference("inputs/phix/PhiX.fa")
	settings := gobwa.GoBwaAllocSettings()
	arena := gobwa.NewArena()
	fastq, _ := fastqreader.OpenFastQ("inputs/zero_length_read_test.fastq.gz")
	improper := -17.0
	inference.SetArgsForTests(inference.LariatArgs{Improper_pair_penalty: &improper})
	bc_reads := make([]fastqreader.FastQRecord, 0, 50000)
	bc_reads, _, _ = fastq.ReadBarcodeSet(&bc_reads, 7)
	barcode_chains, _ := inference.GetChains(ref, settings, bc_reads, arena, 25)
	inference.GetAlignments(ref, settings, barcode_chains, 17.0, arena)
	// if it did not crash, it passes
}
