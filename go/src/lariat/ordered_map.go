// Copyright (c) 2015 10X Genomics, Inc. All rights reserved.

package main

type OrderedMap struct {
	index         map[int]int //index of a key k
	reverse_index []int       //key that lives at an index i
	store         []interface{}
}

func NewOrderedMap() *OrderedMap {
	om := &OrderedMap{
		index:         make(map[int]int),
		reverse_index: []int{},
		store:         make([]interface{}, 0),
	}
	return om
}

func (om *OrderedMap) Get(key int) interface{} {
	toRet, ok := om.index[key]
	if ok {
		return om.store[toRet]
	}
	return nil
}

func (om *OrderedMap) Set(key int, val interface{}) {
	i, ok := om.index[key]
	if ok {
		om.store[i] = val
	} else {
		om.index[key] = len(om.store)
		om.reverse_index = append(om.reverse_index, key)
		om.store = append(om.store, val)
	}
}

func (om *OrderedMap) Delete(key int) {
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

func FixGetForTypeAlignment(val interface{}) *Alignment {
	if val == nil {
		return nil
	}
	return val.(*Alignment)
}

func FixGetForTypeOrderedMap(val interface{}) *OrderedMap {
	if val == nil {
		return nil
	}
	return val.(*OrderedMap)
}

func (om *OrderedMap) Iter() []interface{} {
	return om.store
}

func (om *OrderedMap) IterKeys() []int {
	return om.reverse_index
}

func (om *OrderedMap) Len() int {
	return len(om.reverse_index)
}
