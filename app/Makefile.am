## Process this file with automake to produce Makefile.in

## Created by Anjuta

AM_CPPFLAGS =\
	-DPACKAGE_LOCALE_DIR=\""$(localedir)"\" \
	-DPACKAGE_SRC_DIR=\""$(srcdir)"\" \
	-DPACKAGE_DATA_DIR=\""$(pkgdatadir)"\" \
	$(CGPAPP_CFLAGS)

AM_CFLAGS =\
	-I$(top_builddir)/src\
	 -Wall\
	 -g

bin_PROGRAMS =\
	cgpapp

cgpapp_SOURCES =\
	main.c

cgpapp_LDFLAGS =\
	-lm 

cgpapp_LDADD =\
	$(CGPAPP_LIBS)\
	$(top_builddir)/src/libcgp.la -lm

