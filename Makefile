########################################################################
##
## -----------------------
##   Makefile for KalDet
## -----------------------
##
## (Update Record)
##   2009/11/23  K.Ikematsu   Derived from yokaRawMonNeo/kalman/Makefile
##
## (Description)
##   In order to use this package you should first set some
##   environmental variables. Take a look at README in the
##   top directory.
##
## (Targets)
##   all           : builds libraries
##   Makefiles     : creates Makefiles
##   clean         : deletes *.o ...
##   distclean     : deletes even libraries
##
## $Id: Makefile,v 1.1.1.1 2009-11-24 00:13:59 ikematsu Exp $
########################################################################

SHELL   = `which bash`
MFLAGS  =
CURRDIR = .

#SUBDIRS = kern lctpc/lp1jgem ild/tpc othertpc/mp-tpc
SUBDIRS = kern lctpc/lp1jgem

all:
	@case '${MFLAGS}' in *[ik]*) set +e;; esac; \
	for i in $(SUBDIRS); do \
	(cd $$i; echo ``making'' all ``in $(CURRDIR)/$$i...''; \
	xmkmf -a; \
	$(MAKE) $(MFLAGS) install); \
	done
	@case '${MFLAGS}' in *[ik]*) set +e;; esac; \

Makefiles:
	@case '${MFLAGS}' in *[ik]*) set +e;; esac; \
	for i in $(SUBDIRS); do \
	(cd $$i; echo ``making'' Makefiles ``in $(CURRDIR)/$$i...''; \
	xmkmf); \
	done

clean: Makefiles
	@case '${MFLAGS}' in *[ik]*) set +e;; esac; \
	for i in $(SUBDIRS); do \
	(cd $$i; echo ``making'' clean ``in $(CURRDIR)/$$i...''; \
	$(MAKE) $(MFLAGS) clean); \
	done

distclean: clean
	rm -rf include lib
	@case '${MFLAGS}' in *[ik]*) set +e;; esac; \
	for i in $(SUBDIRS); do \
	(cd $$i; echo ``making'' distclean ``in $(CURRDIR)/$$i...''; \
	$(MAKE) $(MFLAGS) distclean); \
	done; \
