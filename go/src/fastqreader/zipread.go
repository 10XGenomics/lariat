// Copyright (c) 2015 10X Genomics, Inc. All rights reserved.

package fastqreader

import (
	"io"
	"os"
	"os/exec"
)

/*
 * We need a "ZipReader" reader object. A ZipReader contains a source reader
 * but implements a read function that retries if the source reader returned
 * fewer bytes than requested.
 */
type ZipReader struct {
	/* The exact reader that was passed to us */
	Source io.Reader
	Cmd    *exec.Cmd
}

/*
 * Implement a "ReRead" function. It reads from a Reader and it retries
 * if it gets fewer bytes than it wanted.  It will always completely fill
 * the "data" array or return an error.
 */
func (r ZipReader) Read(data []byte) (int, error) {

	/* How much data do we wish to read */
	want := len(data)
	/* Did something go wrong? */
	var err error
	/* Where are we currently reading to in data */
	var offset int
	/* How much data did we just read*/
	var read_len int
	/* Real programmers never use curley braces. They just write all their code as
	 * increment and check conditions in a body-less for loop.
	 */
	for offset = 0; offset < want && err == nil; read_len, err = r.Source.Read(data[offset:]) {
		offset += read_len
	}

	return offset, err
}

func (r ZipReader) Close() {

	r.Source.(io.ReadCloser).Close()
	r.Cmd.Wait()
}

func MakeZipReader(src io.Reader, cmd *exec.Cmd) *ZipReader {
	return &ZipReader{src, cmd}
}

/*
 * Implement an "uncompressing" reader that uses the system's
 * gunzip program to uncompress the input file.  We do this
 * because the system gzip is much much faster than the go library
 */
func FastZipReader(path string) (*ZipReader, error) {
	/* Setup a "Cmd" structure to describe the command we
	 * want to execute
	 * */
	if _, err := os.Stat(path); os.IsNotExist(err) {
		panic("no such file or directory: " + path)
	}
	cmd := exec.Command("gunzip", "-c", path)
	stdout, err := cmd.StdoutPipe()

	if err != nil {
		return nil, err
	}
	/* Start the command */
	err = cmd.Start()

	if err != nil {
		return nil, err
	}

	/* Make a reader that will return data with no short reads from
	 * gzip's stdout */
	return MakeZipReader(stdout, cmd), nil
}
