// Copyright (c) 2015 10X Genomics, Inc. All rights reserved.

/*
THis is a simple example that demonstrates how to use BWA frmo go.
*/

package main

import "bufio"
import "os"
import . "gobwa"
// import "log"
import "fmt"

func main() {
	print("Waking up\n")
	ref := GoBwaLoadReference("/mnt/opt/refdata/fasta/hg19/hg19.fa")
	print("Reference loaded\n")
	settings := GoBwaAllocSettings()
	arena := NewArena()
	stdin := bufio.NewReader(os.Stdin)
	for {
		print("> ")
		read1, _, _ := stdin.ReadLine()
		read1s := string(read1)
		print("> ")
		read2, _, _ := stdin.ReadLine()


					read1 = []byte(read1s)
					fmt.Println(string(read1))
					fmt.Println(string(read2))
					//  chains := GoBwaChain(ref, settings, string(line))
					      //  alignments := GoBwaAlign(ref, settings, string(line), arena)
								 alignments1, alignments2 := GoBwaMemMateSW(ref , settings , &read1 , &read2 , arena , 25)
								//  log.Printf("%v",alignments)
								for i := 0; i < len(alignments1); i++ {
									// log.Printf("%v",alignments[i])
									fmt.Println(alignments1[i].Contig, alignments1[i].Offset,alignments1[i].Score, alignments1[i].Reversed)
								}
								fmt.Println("next")
								for i := 0; i < len(alignments2); i++ {
									// log.Printf("%v",alignments[i])
									fmt.Println(alignments2[i].Contig, alignments2[i].Offset,alignments2[i].Score, alignments2[i].Reversed)
								}


								//  log.Printf("%v",alignments)

		// print("Processing....")


			// fmt.Println("num",len(alignments))
			// for i := 0; i < len(alignments); i++ {
			// 	fmt.Println("chr",alignments[i].Contig,"offset",alignments[i].Offset,"score",alignments[i].Score)
			// 	fmt.Println(string(line2))
			// 		// GoBwaMemMateSmithWaterman(ref, settings, &alignments[i], &line2, arena)
			// }
        // for i := 0; i < len(alignments); i++ {
        //     GoBwaSmithWaterman(ref, settings, string(line), alignments[i])
        // }
	}
}
