
lib_LTLIBRARIES = \
	libcgp.la

libcgp_la_SOURCES = \
	cgp.c \
	cgp.h

libcgp_la_LIBADD = \
	$(CGPAPP_LIBS) \
	-lm

AM_CPPFLAGS = \
	$(CGPAPP_CFLAGS) \
	-DG_LOG_DOMAIN=\"CGPAPP\"

libcgpinclude_HEADERS = \
	cgp.h

EXTRA_DIST = \
	$(libcgpinclude_HEADERS)\
	$(pkgconfig_DATA) 

libcgpincludedir = \
	$(includedir)/libcgp

pkgconfig_DATA = \
	libcgp-1.0.pc

pkgconfigdir = \
	$(libdir)/pkgconfig

## File created by the gnome-build tools

