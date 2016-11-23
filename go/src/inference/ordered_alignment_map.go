// Copyright (c) 2015 10X Genomics, Inc. All rights reserved.

package inference

type OrderedAlignmentMap struct {
	index         map[int]int //index of a key k
	reverse_index []int       //key that lives at an index i
	store         []*Alignment
}

func NewOrderedAlignmentMap() *OrderedAlignmentMap {
	om := &OrderedAlignmentMap{
		index:         make(map[int]int),
		reverse_index: []int{},
		store:         make([]*Alignment, 0),
	}
	return om
}

func (om *OrderedAlignmentMap) Get(key int) *Alignment {
	toRet, ok := om.index[key]
	if ok {
		return om.store[toRet]
	}
	return nil
}

func (om *OrderedAlignmentMap) Set(key int, val *Alignment) {
	i, ok := om.index[key]
	if ok {
		om.store[i] = val
	} else {
		om.index[key] = len(om.store)
		om.reverse_index = append(om.reverse_index, key)
		om.store = append(om.store, val)
	}
}

func (om *OrderedAlignmentMap) Delete(key int) {
	i, ok := om.index[key]
	if ok {
		if len(om.store) > 1 {
			om.store[i] = om.store[len(om.store)-1]
			om.index[om.reverse_index[len(om.store)-1]] = i
			om.reverse_index[i] = om.reverse_index[len(om.reverse_index)-1]
		}
		om.store = om.store[0 : len(om.store)-1]
		om.reverse_index = om.reverse_index[0 : len(om.reverse_index)-1]
		delete(om.index, key)
	}
}

func (om *OrderedAlignmentMap) Iter() []*Alignment {
	return om.store
}

func (om *OrderedAlignmentMap) IterKeys() []int {
	return om.reverse_index
}

func (om *OrderedAlignmentMap) Len() int {
	return len(om.reverse_index)
}
