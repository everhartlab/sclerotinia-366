#!/bin/bash
#
# 2017-11-17
#
# This script is an attempt at addressing some of the issues brought up by
# PeerJ to address their workflow. Specifically we need to:
#
# 1. Separate out the supplementary figures/tables from the manuscript
# 2. Rename figures as Fig1.pdf ... Fig6.pdf for uploading.
#
# This script will perform all steps, including rendering to PDF.
TEXFILE="doc/manuscript/manuscript.tex"

# step 1: figure out where the supplementary information begins and ends.
supplement_line=$(grep -n 'section\*{Supplementary' ${TEXFILE} | awk -F: '{print $1}')
references_line=$(grep -n 'section\*{References' ${TEXFILE} | awk -F: '{print $1}')

echo "Supplementary info starts on line ${supplement_line} and ends on $((references_line - 1))"

# step 2: find the figures
#    - Use () to create an array from the numbers
figures=($(grep -n 'includegraphics' ${TEXFILE} | awk -F: '{print $1}'))
#    - Get six elements of the array, starting at index 0: https://stackoverflow.com/a/1336245/2752888
figures=(${figures[@]:0:6})

echo "Figures are located on lines $(printf '%s ' ${figures[*]})"
echo "I am updating their names now ..."
counter=1;
for i in ${figures[@]};
do
  fig=$(sed "${i}q;d" ${TEXFILE} | sed -r "s/^.+\{(.+)\}/\1/")
  new=$(sed -r "s@^(.+)/([^/]+\.pdf)@\1/Fig${counter}.pdf@" <<< ${fig})
  echo "Figure ${fig} will become ${new}"
  counter=$((counter + 1))
done
