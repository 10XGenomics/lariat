// Copyright (c) 2015 10X Genomics, Inc. All rights reserved.

package test

import "testing"

func Check(t *testing.T, mustbetrue bool, message string) {
	if !mustbetrue {
		t.Error(message)
	}
}
