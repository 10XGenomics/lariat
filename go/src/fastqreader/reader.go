// Copyright (c) 2015 10X Genomics, Inc. All rights reserved.

package fastqreader

import (
	"bufio"
	"io"
	"log"
	"strings"
)

/*
 * This structure represents a single read from a fastq file
 */
type FastQRecord struct {
	Read1     []byte
	ReadQual1 []byte
	Read2     []byte
	ReadQual2 []byte

	Barcode10X     []byte
	Barcode10XQual []byte

	Barcode     []byte
	BarcodeQual []byte

	ReadInfo string
}

/*
 * This struture reprensets a "fastQ" reader that can pull single records
 * as well as sets of records (on the same barcode) from a fastq file
 */
type FastQReader struct {
	Source        *ZipReader
	Buffer        *bufio.Reader
	Line          int
	DefferedError error
	Pending       *FastQRecord
	LastBarcode   *FastQRecord
}

/*
 * A utility function to compare two slices
 */
func SliceCompare(a []byte, b []byte) bool {
	if len(a) != len(b) {
		return false
	}

	for i := 0; i < len(a); i++ {
		if a[i] != b[i] {
			return false
		}
	}
	return true

}

/* Open a new fastQ files */
func OpenFastQ(path string) (*FastQReader, error) {

	var res = new(FastQReader)

	var err error
	res.Source, err = FastZipReader(path)

	if err != nil {
		return nil, err
	}

	res.Buffer = bufio.NewReader(res.Source)
	res.Line = 0
	return res, nil
}

/*
 * Read a single record from a fastQ file
 */
func (fqr *FastQReader) ReadOneLine(result *FastQRecord) error {

	/* Search for the next start-of-record.*/
	for {
		fqr.Line++
		line, err := fqr.Buffer.ReadString(byte('\n'))
		if err != nil {
			return err
		}
		if line[0] == byte('@') {
			/* Found it! */
			//result.ReadInfo = strings.TrimSpace(string(line))
			result.ReadInfo = strings.Fields(string(line[1 : len(line)-1]))[0]
			break
		} else {
			log.Printf("Bad line: %v at %v", string(line), fqr.Line)
		}
	}

	/* Load the 8 sequences for this record */
	var stuff_to_get [8][]byte

	for i := 0; i < 8; i++ {
		var err error
		var line []byte
		line, err = fqr.Buffer.ReadBytes(byte('\n'))
		stuff_to_get[i] = line[0 : len(line)-1]
		if err != nil {
			return err
		}
	}

	/* Assign them to the right fields in the FastQRecord struct */
	result.Read1 = stuff_to_get[0]
	result.ReadQual1 = stuff_to_get[1]
	result.Read2 = stuff_to_get[2]
	result.ReadQual2 = stuff_to_get[3]
	result.Barcode10X = stuff_to_get[4]
	result.Barcode10XQual = stuff_to_get[5]
	result.Barcode = stuff_to_get[6]
	result.BarcodeQual = stuff_to_get[7]

	return nil
}

/*
 * Decide of two reads come from different gems
 */
func DifferentBarcode(a *FastQRecord, b *FastQRecord) bool {
	if SliceCompare(a.Barcode10X, b.Barcode10X) {
		return false
	} else {
		return true
	}
}

/*
 * Reaturn an array of all of the reads from the same GEM.
 * "space" may be null or may be the result of a previous call to this function.
 * If present the array will be destructively re-used
 */
func (fqr *FastQReader) ReadBarcodeSet(space *[]FastQRecord) ([]FastQRecord, error, bool) {
	new_barcode := false
	if fqr.DefferedError != nil {
		return nil, fqr.DefferedError, false
	}

	var record_array []FastQRecord
	if space == nil {
		/* Allocate some space, guessing at most 1 million reads per
		 * barcode. GO will transparently extend this array if needed
		 */
		record_array = make([]FastQRecord, 0, 1000000)
	} else {
		/* Re-use (but truncate) space */
		record_array = (*space)[0:0]
	}

	var index = 0

	/* Is there a pending element from a previous call that needs to be
	 * put in the output?
	 */
	if fqr.Pending != nil {
		record_array = append(record_array, *fqr.Pending)
		fqr.Pending = nil
		index++
	}

	/* Load fastQ records into record_array */
	for ; index < 30000; index++ {
		record_array = append(record_array, FastQRecord{})
		err := fqr.ReadOneLine(&record_array[index])

		if err != nil {
			/* Something went wrong. If we have data, return it and
			 * defer the error to the next invocation. Otherwise,
			 * return the error now.
			 */
			if err != io.EOF {
				log.Printf("Error: %v", err)
			}
			if index == 0 {
				return nil, err, false
			} else {
				fqr.DefferedError = err
				break
			}
		}

		if DifferentBarcode(&record_array[0], &record_array[index]) {
			/* Just transitioned to a new GEM. This record needs to
			 * be defered for next time we're called (since its on the
			 * _new_ gem).
			 */
			fqr.Pending = new(FastQRecord)
			*fqr.Pending = record_array[index]
			new_barcode = true
			break
		} else if fqr.LastBarcode != nil && !DifferentBarcode(&record_array[0], fqr.LastBarcode) && index >= 1000 {
			new_barcode = false
			break
		}

	}
	if len(record_array) > 0 {
		fqr.LastBarcode = &record_array[0]
	}

	//log.Printf("Load %v record %s %s %s %s", index, string(record_array[0].Barcode10X), string(record_array[index].Barcode10X), string(record_array[0].Barcode), string(record_array[index].Barcode))
	/* Truncate the last record of the array. It is either eroneous and ill defined
	 * or it belongs to the next GEM.
	 */

	end := len(record_array)
	if new_barcode {
		end -= 1
	} else if fqr.DefferedError != io.EOF {
		return record_array[0:end], nil, false
	}
	return record_array[0:end], nil, true

}
