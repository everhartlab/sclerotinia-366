.PHONY: all

all: bootstrap \
	results/data-comparison.md \
	results/MLG-distribution.md \
	results/mlg-mcg.md \
	results/RDA-analysis.md \
	results/pop-diff.md

.PHONY: bootstrap

bootstrap: results/bootstrap.txt

results/bootstrap.txt: DESCRIPTION
	R --slave -e "devtools::install()"
	date > results/bootstrap.txt

results/%.md : doc/RMD/%.Rmd
	R --slave -e "ezknitr::ezknit('$^', \
	              fig_dir = './figures/$*/', \
	              out_dir = '$(@D)', \
	              keep_html = FALSE \
	              )"
