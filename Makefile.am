## Process this file with automake to produce Makefile.in
## Created by Anjuta

SUBDIRS = src po app

dist_doc_DATA = \
	README.md \
	LICENSE \
	ChangeLog \
	INSTALL


INTLTOOL_FILES = intltool-extract.in \
	intltool-merge.in \
	intltool-update.in

EXTRA_DIST = \
	$(INTLTOOL_FILES)

DISTCLEANFILES = intltool-extract \
	intltool-merge \
	intltool-update \
	po/.intltool-merge-cache


# Remove doc directory on uninstall
uninstall-local:
	-rm -r $(docdir)
