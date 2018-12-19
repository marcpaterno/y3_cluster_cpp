DATAFILES=$(wildcard ../data/*.out ../data/*.txt)

all: physics_validation.html

physics_validation.md: physics_validation.Rmd $(DATAFILES)
	rscript -e "knitr::knit('$<', '$@')"

physics_validation.html: physics_validation.md
	pandoc --from markdown +RTS -K512M -RTS $< \
		--to html4 \
		--from markdown+autolink_bare_uris+ascii_identifiers+tex_math_single_backslash+smart \
		--output $@ \
		--email-obfuscation none \
		--self-contained \
		--standalone \
		--section-divs \
		--template /Library/Frameworks/R.framework/Versions/3.5/Resources/library/rmarkdown/rmd/h/default.html \
		--no-highlight \
		--variable highlightjs=1 \
		--variable 'theme:bootstrap' \
		--mathjax \
		--variable 'mathjax-url:https://mathjax.rstudio.com/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML'

clean:
	rm -r *.md *.html
	rm -f figure/*
