DATAFILES=$(wildcard ../data/*.out ../data/*.txt)

all: physics_validation.html

physics_validation.html: physics_validation.Rmd $(DATAFILES)
	rscript -e "rmarkdown::render('$<', 'html_document', '$@', clean = FALSE)"

clean:
	rm -f *.md *.html
	rm -f figure/*
