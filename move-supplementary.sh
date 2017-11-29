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
SUPFILE="doc/manuscript/supplementary.tex"
PAPER="doc/manuscript/paper.tex"
MANUSCRIPT_SHELL="doc/manuscript/manuscript_shell.tex"
PAPER_ONLY="doc/manuscript/paper-only.tex"
SUPPLEMENTARY_ONLY="doc/manuscript/supplementary-only.tex"

# step 1: figure out where the supplementary information begins and ends.
supplement_line=$(grep -n 'section\*{Supplementary' ${TEXFILE} | awk -F: '{print $1}')
references_line=$(grep -n 'section\*{References' ${TEXFILE} | awk -F: '{print $1}')
intro_line=$(grep -n 'section\*{Introduction' ${TEXFILE} | awk -F: '{print $1}')
end_line=$(grep -n '\end{document' ${TEXFILE} | awk -F: '{print $1}')

printf "Supplementary info starts on line ${supplement_line} and ends on $((references_line - 1))\n"

# step 2: find the figures
#    - Use () to create an array from the numbers
figures=($(grep -n 'includegraphics' ${TEXFILE} | awk -F: '{print $1}'))
#    - Get six elements of the array, starting at index 0: https://stackoverflow.com/a/1336245/2752888
figures=(${figures[@]:0:6})

printf "\nFigures are located on lines $(printf '%s ' ${figures[*]})\n"
if [ $(ls results/figures/publication | grep -c "Fig1.pdf") -eq 0 ]; then
  printf "I am updating their names now ...\n"
  counter=1;
  for i in ${figures[@]};
  do
    fig=$(sed "${i}q;d" ${TEXFILE} | sed -r "s/^.+\{(.+)\}/\1/" | sed 's@../../@@')
    new=$(sed -r "s@^(.+)/([^/]+\.pdf)@\1/Fig${counter}.pdf@" <<< ${fig})
    printf "\tsed -i \"s@${fig}@${new}@\" ${TEXFILE}\n"
    sed -i "s@${fig}@${new}@" ${TEXFILE}
    printf "\tcp ${fig} ${new}\n"
    cp ${fig} ${new}
    now=$(sed "${i}q;d" ${TEXFILE} | sed -r "s/^.+\{(.+)\}/\1/" | sed 's@../../@@')
    printf "\tFigure ${fig} is now ${now}\n\n"
    counter=$((counter + 1))
  done
else
  printf "Figures are up to date.\n"
fi;

# step 3: split the tex documents into two files to follow https://tex.stackexchange.com/a/87011/77699
printf "\nSplitting the document...\n"

split_sup="sed -n ${supplement_line},$((references_line-1))p ${TEXFILE} > ${SUPFILE}"
printf "\t${split_sup}\n"
eval ${split_sup}

split_body="sed -n ${intro_line},$((supplement_line-1))p ${TEXFILE} > ${PAPER}"
printf "\t${split_body}\n"
eval ${split_body}

split_ref="sed -n ${references_line},$((end_line-1))p ${TEXFILE} >> ${PAPER}"
printf "\t${split_ref}\n"
eval ${split_ref}

printf "\nCreating new document...\n"

split_head="sed -n 1,$((intro_line-1))p ${TEXFILE} > ${MANUSCRIPT_SHELL}"
printf "\t${split_head}\n"
eval ${split_head}

printf "\\\\input{paper}" >> ${MANUSCRIPT_SHELL}
printf "\n\\\\newpage\n" >> ${MANUSCRIPT_SHELL}
printf "\\\\input{supplementary}" >> ${MANUSCRIPT_SHELL}
printf "\n" >> ${MANUSCRIPT_SHELL}

split_foot="printf '\\\\end{document}' >> "
printf "\t${split_foot}${MANUSCRIPT_SHELL}\n"
eval "${split_foot}${MANUSCRIPT_SHELL}"
printf "\n" >> ${MANUSCRIPT_SHELL}


printf "\nCreating manuscript...\n"
cd doc/manuscript
pdflatex manuscript_shell
pdflatex manuscript_shell
pdflatex manuscript_shell
cd -

printf "\ndone\n"
printf "\n\nNow open manuscript_shell.pdf and split it into the supplementary info and the manuscript manually\n"
