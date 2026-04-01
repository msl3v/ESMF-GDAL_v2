# $Id$

ALL:      build_examples
run:      run_examples
run_uni:  run_examples_uni

LOCDIR	  = src/Infrastructure/LineStream/examples


EXAMPLES_BUILD    = $(ESMF_EXDIR)/ESMF_LineStreamEx


EXAMPLES_RUN      = run_ESMF_LineStreamEx

EXAMPLES_RUN_UNI  = run_ESMF_LineStreamEx_uni


include $(ESMF_DIR)/makefile

DIRS	  = 

CLEANDIRS   = 
CLEANFILES  = $(EXAMPLES_BUILD)
CLOBBERDIRS =

#
#  ESMF_LineStream
#
run_ESMF_LineStreamEx:
	$(MAKE) EXNAME=LineStream NP=4 exfrun

run_ESMF_LineStreamEx_uni:
	$(MAKE) EXNAME=LineStream NP=1 exfrun





