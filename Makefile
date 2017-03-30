all: results/data-comparison.Rmd

.PHONY: all

results/%.md : doc/RMD/%.Rmd
R --slave -e "ezknitr::ezknit('$^', \
              fig_dir = './figures/$*/', \
              out_dir = '$(@D)', \
              keep_html = FALSE \
              )"
