SHELL := /bin/sh
CC := gcc
CFLAGS := -Wextra -Wpedantic -Wall -std=c17

INSTALL := install
INSTALL_PROGRAM := $(INSTALL)
INSTALL_DATA := $(INSTALL) -m 644

prefix := /usr/local
exec_prefix := $(prefix)
bindir := $(exec_prefix)/bin

OS := $(shell uname -s)
ifeq ($(OS),Darwin)
	CC := clang
	CPPFLAGS += -I/opt/homebrew/include
	LDFLAGS += -L/opt/homebrew/lib
endif

.SUFFIXES : .c .o

.PHONY: all
all: release

.PHONY: filtergt
filtergt: filtergt.c
	$(CC) $(CPPFLAGS) $(CFLAGS) $(LDFLAGS) $< -lhts $(LDLIBS) -o $@

.PHONY: release
release: CFLAGS := -O2 $(CFLAGS)
release: filtergt

.PHONY: debug
debug: CFLAGS := -g $(CFLAGS)
debug: filtergt

.PHONY: install
install: filtergt
	$(INSTALL_PROGRAM) $< $(DESTDIR)$(bindir)/$<

.PHONY: uninstall
uninstall: $(DESTDIR)$(bindir)/filtergt
	rm -f $(DESTDIR)$(bindir)/$<
