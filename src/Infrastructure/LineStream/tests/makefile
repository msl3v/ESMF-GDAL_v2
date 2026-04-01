# $Id$ 

ALL:      build_unit_tests
run:      run_unit_tests
run_uni:  run_unit_tests_uni

LOCDIR	  = src/Infrastructure/LineStream/tests

.NOTPARALLEL:
TESTS_BUILD   = \
	$(ESMF_TESTDIR)/ESMF_LineStreamUTest \
	$(ESMF_TESTDIR)/ESMC_LineStreamUTest


TESTS_RUN     = \
	RUN_ESMF_LineStreamUTest \
	RUN_ESMC_LineStreamUTest 

TESTS_RUN_UNI = \
	RUN_ESMF_LineStreamUTestUNI \
	RUN_ESMC_LineStreamUTestUNI 

include ${ESMF_DIR}/makefile

CLEANDIRS   = 
CLEANFILES  = $(TESTS_BUILD)
CLOBBERDIRS =

DIRS        = 

#
# unit test targets
#

# --- LineStream Creation

RUN_ESMF_LineStreamUTest:
	$(MAKE) TNAME=LineStream NP=4 ftest

RUN_ESMF_LineStreamUTestUNI:
	$(MAKE) TNAME=LineStream NP=1 ftest

RUN_ESMC_LineStreamUTest:
	$(MAKE) TNAME=LineStream NP=4 ctest

RUN_ESMC_LineStreamUTestUNI:
	$(MAKE) TNAME=LineStream NP=1 ctest
