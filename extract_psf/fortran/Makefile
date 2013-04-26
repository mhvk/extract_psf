all: extrpsf.so

#extract_psf.f95: $(MARTEN_GENERAL)/reduce/extract_psf_new.f95
#	cp $< $@
#	echo "subroutine writeout(out)\n  character :: out*(*)\n  print*,out\nend" >>$@

extrpsf.so: extract_psf.f95
	f2py -m $(subst .so,,$@) only: extract_psf psfmany : -c $<
