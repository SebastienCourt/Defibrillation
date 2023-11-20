# SUBDIRS = 

check_PROGRAMS = SFB_maxmax SFB_maxmax2

CLEANFILES = 

SFB_maxmax_SOURCES = SFB_maxmax.cc SFB_maxmax.h SFB_maxmax_routines.cc SFB_maxmax_routines.h SFB_maxmax_solving.cc SFB_maxmax_solving.h
SFB_maxmax2_SOURCES = SFB_maxmax2.cc SFB_maxmax2.h SFB_maxmax_routines2.cc SFB_maxmax_routines2.h SFB_maxmax_solving2.cc SFB_maxmax_solving2.h

MUMPS_LIBS = @MUMPS_LIBS@
# MUMPS_CFLAGS = @MUMPS_CFLAGS@
AM_CPPFLAGS = -I$(top_srcdir)/src -I../../src # $(MUMPS_CFLAGS)
LDADD    = ../../src/libgetfem.la -lm @SUPLDFLAGS@



#TESTS = shallow_water.pl

EXTRA_DIST = \
	SFB_maxmax.param	\
	data_elasticity.param		

LOG_COMPILER = perl
