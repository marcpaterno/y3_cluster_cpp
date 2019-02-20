DATAFILES=$(wildcard ../data/*.out ../data/*.txt)

all: physics_validation.html

physics_validation.html: physics_validation.Rmd $(DATAFILES)
	rscript -e "rmarkdown::render('$<', 'html_document', '$@', clean = FALSE)"

physics_validation.pdf: physics_validation.Rmd $(DATAFILES)
	rscript -e "rmarkdown::render('$<', 'pdf_document', '$@', clean = FALSE)"

physics_validation.docx: physics_validation.Rmd $(DATAFILES)
	rscript -e "rmarkdown::render('$<', 'word_document', '$@', clean = FALSE)"

pdf: physics_validation.pdf
word: physics_validation.docx

clean:
	rm -f *.md *.html *.pdf *.docx
	rm -f figure/*
	rm -rf physics_validation_files/

