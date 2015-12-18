// Copyright (c) 2015 10X Genomics, Inc. All rights reserved.

package test

import "fastqreader"
import "testing"

func TestReader(t *testing.T) {

	a, b := fastqreader.OpenFastQ("./inputs/1.fq")
	Check(t, a != nil, "1")
	Check(t, b == nil, "2")

	var ff fastqreader.FastQRecord

	a.ReadOneLine(&ff)
	a.ReadOneLine(&ff)
	a.ReadOneLine(&ff)
	Check(t, string(ff.Read1) == "CACCGCCCTAGCCAGGAGAGAAGCACTTCTTACCTGGGTTTCTTAGAGGCTTTGGCTGGCAATATTGTCAGCACCAGAGAGGACTTCTCGATGGCTGA", "a")
	Check(t, string(ff.ReadQual1) == "BBBFFFFFFFFFFIIIIIFFIIIIIIIIFIIIIIFIFIFFIIFIIIIIIIIIIIIIIIFFFFFFFFFFFFFFFFFFFBFFFFFFFFFFFFFFFFFFFF", "b")
	Check(t, string(ff.Read2) == "GTGGTAGTCTCCTGTTCAGCCATCGAGAAGTCCTCTCTGGTGCTGACAATATTGCCAGCCAAAGCCTCTAAGAAACCCAGGTAAGAAGTGCTTCTCTC", "c")
	Check(t, string(ff.ReadQual2) == "BBBFFFFFFFFFFIIFIIFIIIIIIIIIIIIIIIIFIIIFFFIIIIIIIIIIIIIIIIIIIIFIIIIIIIIIIFFFFFFFFBBFFFFFBBFFFFFFFF", "d")
	Check(t, string(ff.Barcode10X) == "AAACAGAGAAAGAT", "e")
	Check(t, string(ff.Barcode10XQual) == "BBBFFFFFFFFFFI", "f")
	Check(t, string(ff.Barcode) == "CCGAACGC", "g")
	Check(t, string(ff.BarcodeQual) == "BBBFFFFF", "h")
	Check(t, ff.ReadInfo == "@HWI-D00684:80:HFCKCADXX:2:2113:9410:56703 1:N:0:\n", "i")
}

func TestReader2(t *testing.T) {
	a, _ := fastqreader.OpenFastQ("./inputs/1.fq")

	set1, _ := a.ReadBarcodeSet(nil)

	Check(t, set1[0].ReadInfo == "@HWI-D00684:80:HFCKCADXX:2:2113:17628:14813 1:N:0:\n", "a")
	Check(t, string(set1[1].Read1) == "AGCTGCTGCTCTCTCCATGTTTTTCCTGCACTCCTTGCAGGGACCTGAATAGCATGAACTGACTTTTCCTTGACGTAGTTGCTTCGTAGGATACTTCT", "b")

	set2, _ := a.ReadBarcodeSet(&set1)

	Check(t, set2[0].ReadInfo == "@HWI-D00684:80:HFCKCADXX:2:2112:14227:100270 1:N:0:\n", "c")

	Check(t, string(set2[1].Read1) == "CGCGGGCAGCAGCCATGGGATGCAGGACCTGCAGTCCACACATGTCACATGAATCTCCATGGAGAGGCACACAGTTCTCCCCATCTCAGCACTCTCTC", "d")
}

func TestZip(t *testing.T) {

	zr, _ := fastqreader.FastZipReader("./inputs/hello.txt.gz")

	x := make([]byte, 14)

	zr.Read(x)

	Check(t, string(x) == "Hello World!\n\000", "1")

	zr.Close()
}
