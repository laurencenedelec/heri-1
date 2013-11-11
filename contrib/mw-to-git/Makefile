#
# Copyright (C) 2013
#     Matthieu Moy <Matthieu.Moy@imag.fr>
#
# To install, run Git's toplevel 'make install' then run:
#
#   make install

GIT_MEDIAWIKI_PM=Git/Mediawiki.pm
SCRIPT_PERL=git-remote-mediawiki.perl
GIT_ROOT_DIR=../..
HERE=contrib/mw-to-git/

SCRIPT_PERL_FULL=$(patsubst %,$(HERE)/%,$(SCRIPT_PERL))
INSTLIBDIR=$(shell $(MAKE) -C $(GIT_ROOT_DIR)/perl \
                -s --no-print-directory instlibdir)

all: build

install_pm:
	install $(GIT_MEDIAWIKI_PM) $(INSTLIBDIR)/$(GIT_MEDIAWIKI_PM)

build:
	$(MAKE) -C $(GIT_ROOT_DIR) SCRIPT_PERL=$(SCRIPT_PERL_FULL) \
                build-perl-script

install: install_pm
	$(MAKE) -C $(GIT_ROOT_DIR) SCRIPT_PERL=$(SCRIPT_PERL_FULL) \
                install-perl-script

clean:
	$(MAKE) -C $(GIT_ROOT_DIR) SCRIPT_PERL=$(SCRIPT_PERL_FULL) \
                clean-perl-script

perlcritic:
	perlcritic -2 *.perl
