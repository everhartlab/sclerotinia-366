#!/bin/bash

# Test if we have all the results. Note, this only tests if the exist. It does
# not test if they have the same results as the previous files. Some of these
# files will change due to randomness of the analyses, but the overall results
# should not change. 
resfiles=$(grep -E results[^%]+md /analysis/Makefile \
| sed -r "s_^.{14}_/analysis/_" \
| sed 's_\\__')

resfiles+=" /analysis/doc/manuscript/manuscript.pdf"
resfiles+=" /analysis/doc/manuscript/manuscript.tex"

for i in ${resfiles}; 
do 
	md5sum $i
done