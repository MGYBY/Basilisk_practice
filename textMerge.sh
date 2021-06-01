#!/bin/bash
for t in $(seq 0 2 60) 
do
	f1="out-h-${t}"
	f2="out-u-${t}"
	paste $f1 $f2 > "out-comb-${t}" 
done 
