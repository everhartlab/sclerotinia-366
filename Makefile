.PHONY: all

all: bootstrap \
	results/MLG-distribution.md \
	results/mlg-mcg.md \
	results/RDA-analysis.md \
	results/pop-diff.md

.PHONY: bootstrap

bootstrap: results/bootstrap.txt

results/bootstrap.txt: DESCRIPTION
	R --slave -e "devtools::install()"
	date > results/bootstrap.txt

results/data-comparison.md : doc/RMD/data-comparison.Rmd
	R --slave -e "ezknitr::ezknit('$<', \
	              fig_dir = './figures/$*/', \
	              out_dir = '$(@D)', \
	              keep_html = FALSE \
	              )"

data/sclerotinia_16_loci.rda : results/data-comparison.md

results/%.md : doc/RMD/%.Rmd data/sclerotinia_16_loci.rda
	R --slave -e "ezknitr::ezknit('$<', \
	              fig_dir = './figures/$*/', \
	              out_dir = '$(@D)', \
	              keep_html = FALSE \
	              )"

