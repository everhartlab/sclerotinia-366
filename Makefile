.PHONY: all

all: results/data-comparison.md \
	results/MLG-distribution.md

results/%.md : doc/RMD/%.Rmd
	R --slave -e "ezknitr::ezknit('$^', \
	              fig_dir = './figures/$*/', \
	              out_dir = '$(@D)', \
	              keep_html = FALSE \
	              )"
