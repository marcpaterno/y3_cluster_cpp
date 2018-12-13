DATAFILES=$(wildcard ../data/*.out ../data/*.txt)

all: physics_validation.html

physics_validation.md: physics_validation.Rmd $(DATAFILES)
	rscript -e "knitr::knit('$<', '$@')"

physics_validation.html: physics_validation.md
	pandoc --from markdown --to html -o $@ $<

clean:
	rm -r *.md *.html
	rm -f figure/*
