// Copyright (c) 2015 10X Genomics, Inc. All rights reserved.

package test

import "fastqreader"
import "testing"

func TestReader(t *testing.T) {

	a, b := fastqreader.OpenFastQ("./inputs/1.fq")
	Check(t, a != nil, "1")
	Check(t, b == nil, "2")

	var ff fastqreader.FastQRecord

	a.ReadOneLine(&ff, 2)
	a.ReadOneLine(&ff, 2)
	a.ReadOneLine(&ff, 2)
	Check(t, string(ff.Read1) == "CCGCCCTAGCCAGGAGAGAAGCACTTCTTACCTGGGTTTCTTAGAGGCTTTGGCTGGCAATATTGTCAGCACCAGAGAGGACTTCTCGATGGCTGA", "a")
	Check(t, string(ff.ReadQual1) == "BFFFFFFFFFFIIIIIFFIIIIIIIIFIIIIIFIFIFFIIFIIIIIIIIIIIIIIIFFFFFFFFFFFFFFFFFFFBFFFFFFFFFFFFFFFFFFFF", "b")
	Check(t, string(ff.Read2) == "GTGGTAGTCTCCTGTTCAGCCATCGAGAAGTCCTCTCTGGTGCTGACAATATTGCCAGCCAAAGCCTCTAAGAAACCCAGGTAAGAAGTGCTTCTCTC", "c")
	Check(t, string(ff.ReadQual2) == "BBBFFFFFFFFFFIIFIIFIIIIIIIIIIIIIIIIFIIIFFFIIIIIIIIIIIIIIIIIIIIFIIIIIIIIIIFFFFFFFFBBFFFFFBBFFFFFFFF", "d")
	Check(t, string(ff.Barcode10X) == "AAACAGAGAAAGAT", "e")
	Check(t, string(ff.Barcode10XQual) == "BBBFFFFFFFFFFI", "f")
	Check(t, string(ff.Barcode) == "CCGAACGC", "g")
	Check(t, string(ff.BarcodeQual) == "BBBFFFFF", "h")
	Check(t, ff.ReadInfo == "HWI-D00684:80:HFCKCADXX:2:2113:9410:56703", "i")
}

func TestReader2(t *testing.T) {
	a, _ := fastqreader.OpenFastQ("./inputs/1.fq")

	set1, _, _ := a.ReadBarcodeSet(nil, 2)
	Check(t, set1[0].ReadInfo == "HWI-D00684:80:HFCKCADXX:2:2113:17628:14813", "a")
	Check(t, string(set1[1].Read1) == "CTGCTGCTCTCTCCATGTTTTTCCTGCACTCCTTGCAGGGACCTGAATAGCATGAACTGACTTTTCCTTGACGTAGTTGCTTCGTAGGATACTTCT", "b")

	set2, _, _ := a.ReadBarcodeSet(&set1, 2)
	Check(t, set2[0].ReadInfo == "HWI-D00684:80:HFCKCADXX:2:2112:14227:100270", "c")

	Check(t, string(set2[1].Read1) == "CGGGCAGCAGCCATGGGATGCAGGACCTGCAGTCCACACATGTCACATGAATCTCCATGGAGAGGCACACAGTTCTCCCCATCTCAGCACTCTCTC", "d")
}

func TestShortR1(t *testing.T) {
	a, _ := fastqreader.OpenFastQ("./inputs/zero_length_read_test.fastq.gz")

	set1, _, _ := a.ReadBarcodeSet(nil, 7)
	Check(t, len(set1) > 0, "got not reads")
}

func TestZip(t *testing.T) {

	zr, _ := fastqreader.FastZipReader("./inputs/hello.txt.gz")

	x := make([]byte, 14)

	zr.Read(x)

	Check(t, string(x) == "Hello World!\n\000", "1")

	zr.Close()
}
