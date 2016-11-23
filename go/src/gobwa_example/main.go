// Copyright (c) 2015 10X Genomics, Inc. All rights reserved.

/*
This is a simple example that demonstrates how to use BWA from go.
*/

package main

import "bufio"
import "os"
import . "gobwa"
import "log"
import "fmt"
import "strings"

func main() {
	print("Waking up\n")
	ref := GoBwaLoadReference("/mnt/opt/refdata_new/hg19-2.0.0/fasta/genome.fa")
	print("Reference loaded\n")
	settings := GoBwaAllocSettings()
	//ref.GetSeq("chr2",50000,50100,false)
	arena := NewArena()
	stdin := bufio.NewReader(os.Stdin)
	max_alignments := 5
	for {
		print("> ")
		read1, _, _ := stdin.ReadLine()
		print("> ")
		read2, _, _ := stdin.ReadLine()

		read1 = []byte(strings.TrimSpace(string(read1)))
		read2 = []byte(strings.TrimSpace(string(read2)))

		if len(read1) == 0 || len(read2) == 0 {
			fmt.Println("Reads 1 and 2 must not be empty")
			continue
		}

		// fmt.Println(string(read1))
		// fmt.Println(string(read2))
		// chains := GoBwaChain(ref, settings, string(line))
		// alignments := GoBwaAlign(ref, settings, string(line), arena)
		alignments1, alignments2 := GoBwaMemMateSW(ref , settings , &read1 , &read2 , arena , 25)
		// log.Printf("%v",alignments)
		fmt.Println("Read 1")
		for i := 0; i < len(alignments1) && i < max_alignments; i++ {
			// log.Printf("%v",alignments[i])
			fmt.Println(alignments1[i].Contig, alignments1[i].Offset,alignments1[i].Score, alignments1[i].Reversed)
		}
		fmt.Println("Read 2")
		for i := 0; i < len(alignments2) && i < max_alignments; i++ {
			// log.Printf("%v",alignments[i])
			fmt.Println(alignments2[i].Contig, alignments2[i].Offset,alignments2[i].Score, alignments2[i].Reversed)
		}
		// log.Printf("%v",alignments)
		// print("Processing....")
		// fmt.Println("num",len(alignments))
		// for i := 0; i < len(alignments); i++ {
		// fmt.Println("chr",alignments[i].Contig,"offset",alignments[i].Offset,"score",alignments[i].Score)
		// fmt.Println(string(line2))
		// GoBwaMemMateSmithWaterman(ref, settings, &alignments[i], &line2, arena)
		// for i := 0; i < len(alignments); i++ {
		//     GoBwaSmithWaterman(ref, settings, string(line), alignments[i])
		// }
	}
}
