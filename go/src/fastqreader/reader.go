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
    TrimBases []byte
    TrimQuals []byte

	Barcode10X     []byte
	Barcode10XQual []byte
    RawBarcode10X  []byte

	Barcode     []byte
	BarcodeQual []byte

	ReadInfo string
	ReadGroupId string
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
	LastBarcode   []byte
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

func Min(x, y int) int {
    if x < y {
        return x
    }
    return y
}


/*
 * Read a single record from a fastQ file
 */
func (fqr *FastQReader) ReadOneLine(result *FastQRecord, trim int) error {

	/* Search for the next start-of-record.*/
	for {
		fqr.Line++
		line, err := fqr.Buffer.ReadString(byte('\n'))
		if err != nil {
			return err
		}
		if line[0] == byte('@') {
			/* Found it! */
			fields := strings.Fields(string(line[1 : len(line)-1]))
			result.ReadInfo = fields[0]
			if len(fields) < 2 {
				result.ReadGroupId = "" // no RGID found
			} else {
				result.ReadGroupId = fields[len(fields)-1]
			}
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
    var to_trim = Min(len(stuff_to_get[0]), trim)
    r1 := stuff_to_get[0][to_trim:]
    rq := stuff_to_get[1][to_trim:]
    tb := stuff_to_get[0][0:to_trim]
    tq := stuff_to_get[1][0:to_trim]
	result.Read1 = r1
	result.ReadQual1 = rq
    result.TrimBases = tb
    result.TrimQuals = tq
	result.Read2 = stuff_to_get[2]
	result.ReadQual2 = stuff_to_get[3]

    barcodes := strings.Split(string(stuff_to_get[4]), ",")
    result.Barcode10X = []byte(barcodes[0])
    result.RawBarcode10X = []byte(barcodes[len(barcodes)-1])
	result.Barcode10XQual = stuff_to_get[5]
	result.Barcode = stuff_to_get[6]
	result.BarcodeQual = stuff_to_get[7]

	return nil
}

/*
 * Decide of two reads come from different gems
 */
func DifferentBarcode(a []byte, b []byte) bool {
	if SliceCompare(a, b) {
		return false
	} else {
		return true
	}
}

func NotWhitelist(a *FastQRecord) bool {
    for i := 0; i < len(a.Barcode10X); i++ {
        if a.Barcode10X[i] == '-' {
            return false
        }
    }
    return true
}

/*
 * Reaturn an array of all of the reads from the same GEM.
 * "space" may be null or may be the result of a previous call to this function.
 * If present the array will be destructively re-used
 */
func (fqr *FastQReader) ReadBarcodeSet(space *[]FastQRecord, trim int) ([]FastQRecord, error, bool) {
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
		err := fqr.ReadOneLine(&record_array[index], trim)

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

		if DifferentBarcode(record_array[0].Barcode10X, record_array[index].Barcode10X) || (NotWhitelist(&record_array[0]) && index >= 200) {
			/* Just transitioned to a new GEM. This record needs to
			 * be defered for next time we're called (since its on the
			 * _new_ gem).
			 */
			fqr.Pending = new(FastQRecord)
			*fqr.Pending = record_array[index]
			new_barcode = true
			break
		} else if fqr.LastBarcode != nil && !DifferentBarcode(record_array[0].Barcode10X, fqr.LastBarcode) && index >= 200 {
			new_barcode = false
			log.Printf("abnormal break: %s", string(record_array[0].Barcode10X))
			break
		}

	}
	if len(record_array) > 0 {
		tmp := make([]byte, len(record_array[0].Barcode10X))
		copy(tmp, record_array[0].Barcode10X)
		fqr.LastBarcode = tmp
	}
	//log.Printf("Load %v record %s %s %s %s", index, string(record_array[0].Barcode10X), string(record_array[index].Barcode10X), string(record_array[0].Barcode), string(record_array[index].Barcode))
	/* Truncate the last record of the array. It is either eroneous and ill defined
	 * or it belongs to the next GEM.
	 */

	end := len(record_array)
	if new_barcode || fqr.DefferedError == io.EOF {
		end -= 1
	} else if fqr.DefferedError != io.EOF {
		return record_array[0:end], nil, false
	}
	return record_array[0:end], nil, true

}
