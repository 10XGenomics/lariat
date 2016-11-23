// Copyright (c) 2015 10X Genomics, Inc. All rights reserved.

package test

import (
	"gobwa"
	"log"
	"testing"
)

func TestGobwa1(t *testing.T) {

    arena := gobwa.NewArena()
	ref := gobwa.GoBwaLoadReference("inputs/phix/PhiX.fa")
	Check(t, ref != nil, "a")
	settings := gobwa.GoBwaAllocSettings()

	algns := gobwa.GoBwaAlign(ref, settings, "TCAAAAACTGACGCGTTGGATGAGGAGAAGTGGCTTAATATGCTTGGCACGTTCGTCAAGGACTGGTTTA", arena)

	log.Printf("FINAL %v", algns)

	Check(t, algns[0].Offset == 210, "c")
	Check(t, algns[0].Contig == "PhiX", "d")
	algns = gobwa.GoBwaAlign(ref, settings, "TATGACCAGTGTTTCCAGTCCGTTCAGTTGTTGCAGTGGAATAGTCAGGTTAAATTTAATGTGACCGCTT", arena)
	Check(t, len(algns) == 1, "e")
}
