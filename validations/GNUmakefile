DATAFILES=../release-build/data/omega_z_sdss.out \
      ../data/masses_test.out \
		  ../release-build/data/dvdodz_test.out \
		  ../release-build/data/mor_tt_test.out \
		  ../release-build/data/lc_lt_test.out \
		  ../release-build/data/lo_lc_test.out \
		  ../data/integrated_bin_test.out \
		  ../data/integrated_bin_y1_test.out \
		  ../deltasigma/test/xi_nfwonly.out ../deltasigma/test/xi_mm_2halo.out ../deltasigma/test/xi_mm_2halo.out ../deltasigma/test/Sigma_nfwonly.out ../deltasigma/test/Sigma_mm_2halo.out ../deltasigma/test/bias.out \
		  ../release-build/data/Rmis_gammadist_test.out \
		  ../release-build/data/Sigma_crit_test.out

all: html
data: $(DATAFILES)
html: physics_validation.html
pdf: physics_validation.pdf
word: physics_validation.docx

physics_validation.html: physics_validation.Rmd $(DATAFILES)
	Rscript -e "rmarkdown::render('$<', 'html_document', '$@', clean = FALSE)"

physics_validation.pdf: physics_validation.Rmd $(DATAFILES)
	Rscript -e "rmarkdown::render('$<', 'pdf_document', '$@', clean = FALSE)"

physics_validation.docx: physics_validation.Rmd $(DATAFILES)
	Rscript -e "rmarkdown::render('$<', 'word_document', '$@', clean = FALSE)"

clean:
	rm -f *.md *.html *.pdf *.docx
	rm -f figure/*
	rm -rf physics_validation_files/
	rm -f ../data/*.out
	rm -rf ../release-build/data/*.out
	rm -rf ../masses_test_output
	rm -rf ../sdss_analysis_output
	rm -rf ../y1_analysis_output
	rm -rf ../OUTPUT
	rm -rf physics_validation_cache/

../release-build/data/omega_z_sdss.out: ../release-build/test/omega_z_sdss_test
	cd ../release-build/test && ./omega_z_sdss_test

# mfp: It appears to me that cosmosis_tests/mf_tinker_test.py requires the file
# ../release-build/test/test_hmf_z0_z03.dat. The file ../release-build/test/test_hmf_z0_z03.dat has been
# checked into the git repository, so for now we only have our target depend on
# that file.
../data/masses_test.out: ../cosmosis_tests/mf_tinker_test.py ../cosmosis_tests/masses_test.ini ../test/test_hmf_z0_z03.dat ../cosmosis_tests/values/values_comp.ini
	@echo "About to run cosmosis to generate $@"
	cd .. && cosmosis cosmosis_tests/masses_test.ini > /tmp/masses_test.out

../release-build/data/dvdodz_test.out: ../release-build/test/dv_do_dz_t_test
	cd ../release-build/test && ./dv_do_dz_t_test

../release-build/data/mor_tt_test.out: ../release-build/test/mor_t2_test
	cd ../release-build/test && ./mor_t2_test

../release-build/data/lc_lt_test.out: ../release-build/test/lc_lt_t_test
	cd ../release-build/test && ./lc_lt_t_test

../release-build/data/lo_lc_test.out: ../release-build/test/lo_lc_t_test
	cd ../release-build/test/&& ./lo_lc_t_test

../data/integrated_bin_test.out: ../cosmosis_tests/integrated_bin_test.py ../cosmosis_tests/sdss_analysis.ini ../cosmosis_tests/values/values_comp.ini
	@echo "About to run cosmosis to generate $@"
	cd .. && cosmosis cosmosis_tests/sdss_analysis.ini > /tmp/sdss_analysis.out

../data/integrated_bin_y1_test.out: ../cosmosis_tests/integrated_bin_y1_test.py ../cosmosis_tests/y1_analysis.ini ../cosmosis_tests/values/values_comp_y1.ini
	@echo "About to run cosmosis to generate $@"
	cd .. && cosmosis cosmosis_tests/y1_analysis.ini > /dev/null

# The one module ../deltasigma/test/deltasigma_test.py, run by
# ../deltasigma/deltasigma_test.ini, generates several output files. This makes
# the 'make' rule for this target extra complicated, since we don't want to
# re-run the same CosmoSIS workflow several times to generate the same set of
# files. Breaking up the giant one-liner defining the target and dependencies
# seems fragile, so if you try to do so, please test the result carefully!
../deltasigma/test/xi_nfwonly.out ../deltasigma/test/xi_mm_2halo.out ../deltasigma/test/Sigma_nfwonly.out ../deltasigma/test/Sigma_mm_2halo.out ../deltasigma/test/bias.out: ../cosmosis_tests/deltasigma_test.py ../deltasigma/deltasigma_test.ini ../src/modules/deltasigma/deltasigma_values.ini
	@echo "About to run cosmosis to generate $@"
	cd .. && cosmosis deltasigma/deltasigma_test.ini > /dev/null

../release-build/data/Rmis_gammadist_test.out: ../release-build/test/roffset_t_test
	cd ../release-build/test && ./roffset_t_test


../release-build/data/Sigma_crit_test.out: ../release-build/test/sigma_crit_inv_test
	cd ../release-build/test && ./sigma_crit_inv_test
    
