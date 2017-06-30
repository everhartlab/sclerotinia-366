PARSE_DATA := results/data-comparison.md
THE_DATA   := data/clean_data.csv \
              data/raw_data.csv \
              data/sclerotinia_16_loci.rda
ANALYSES   := results/table-1.md \
              results/MCG-virulence.md \
              results/locus-stats.md \
              results/MLG-distribution.md \
              results/mlg-mcg.md \
              results/RDA-analysis.md \
              results/pop-diff.md
MANUSCRIPT := doc/manuscript/manuscript.pdf
COMPONENTS := doc/manuscript/abstract.md \
              doc/manuscript/ssc_bibliography.bib \
              doc/manuscript/wlpeerj.cls


# TARGETS
# ---------------------------------------------------------
.PHONY: all
all: $(ANALYSES) $(MANUSCRIPT)
# In reality $(ANALYSES) -> shared_data -> $(PARSE_DATA) -> bootstrap

# Bootstrap the  data by installing the dependencies
# This is a specific target for the convenience of typing make bootstrap
.PHONY: bootstrap
bootstrap: results/bootstrap.txt

# Create the shared data set
$(THE_DATA) : bootstrap $(PARSE_DATA)

# All the analyses (defined above) depend on the shared data set
$(ANALYSES): $(THE_DATA)

# RECIPES
# ---------------------------------------------------------
results/bootstrap.txt: DESCRIPTION
	R --slave -e "devtools::install()"
	date > results/bootstrap.txt

results/%.md : doc/RMD/%.Rmd
	R --slave -e "ezknitr::ezknit('$<', \
	              fig_dir = './figures/$*/', \
	              out_dir = '$(@D)', \
	              keep_html = FALSE \
	              )"

doc/manuscript/%.pdf : doc/manuscript/%.Rmd $(COMPONENTS) $(ANALYSES)
	R --slave -e "rmarkdown::render('$<')"
