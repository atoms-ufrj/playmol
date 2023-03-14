# Define DEBUG or FAST mode:
#   Build the "fast" version with: `make` or `make DEBUG=0`
#   Build the "debug" version with: `make DEBUG=1`
DEBUG ?= 0

# Installation prefix:
PREFIX ?= /usr/local

# Compilers and their basic options:
FORT ?= gfortran
BASIC_F_OPTS = -march=native -cpp -fmax-errors=1 -Wunused

# Option FAST (default):
FAST_F_OPTS = -O3
#removed o -Ofast flag for suspicious of introducing bugs in the lmp file write step
#https://wiki.gentoo.org/wiki/GCC_optimization
#https://gcc.gnu.org/onlinedocs/gcc-4.6.2/gcc/Optimize-Options.html
#https://stackoverflow.com/questions/14492436/g-optimization-beyond-o3-ofast

# Option DEBUG:
DEBUG_F_OPTS = -Wall -Wno-maybe-uninitialized

# Checks chosen option:
ifeq ($(DEBUG), 1)
  FOPTS = $(BASIC_F_OPTS) $(DEBUG_F_OPTS)
else
  FOPTS = $(BASIC_F_OPTS) $(FAST_F_OPTS)
endif

SRCDIR  = ./src
OBJDIR  = $(SRCDIR)/obj
BINDIR  = ./bin
DOCDIR  = ./docs
GTKDIR  = ./highlight

PACKMOL = ./lib

src = $(addprefix $(SRCDIR)/, $(addsuffix .f90, $(1)))
obj = $(addprefix $(OBJDIR)/, $(addsuffix .o, $(1)))
SRC  = $(call src, mBox mFix mMixingRule mMolecule mParser mString mMath mCodeFlow mGlobal \
                   mPackmol mPlaymol mStruc playmol)
AUX  = $(call src, $(addprefix write_, pdb emdee lammps openmm lammpstrj summary internals)) \
       $(SRCDIR)/elements.inc
OBJ  = $(patsubst $(SRCDIR)/%.f90,$(OBJDIR)/%.o,$(SRC))

all: $(BINDIR)/playmol $(BINDIR)/playmoltools $(BINDIR)/packmol

.PHONY: install clean clean-all doc

doc:
	make -C $(DOCDIR)

clean:
	rm -rf $(OBJDIR) $(BINDIR)

clean-all:
	rm -rf $(OBJDIR) $(BINDIR)
	cd $(PACKMOL) && make clean-all

install:
	cp -f $(BINDIR)/* $(PREFIX)/bin
	bash $(GTKDIR)/install.sh

$(BINDIR)/playmol: $(OBJ) $(PACKMOL)/libpackmol.a
	mkdir -p $(BINDIR)
	$(FORT) $(FOPTS) -J$(OBJDIR) -o $@ $^

$(BINDIR)/playmoltools: $(SRCDIR)/playmoltools.py
	mkdir -p $(BINDIR)
	cp -f $< $@

$(BINDIR)/packmol: $(PACKMOL)/packmol/packmol
	mkdir -p $(BINDIR)
	cp -f $< $@

$(PACKMOL)/libpackmol.a:
	cd $(PACKMOL) && make

$(OBJDIR)/playmol.o: $(SRCDIR)/playmol.f90 $(OBJDIR)/mPlaymol.o
	$(FORT) $(FOPTS) -J$(OBJDIR) -c -o $@ $<

$(OBJDIR)/mPlaymol.o: $(SRCDIR)/mPlaymol.f90 $(AUX) \
                      $(call obj,mCodeFlow mPackmol mFix mMixingRule mBox mMolecule)
	$(FORT) $(FOPTS) -J$(OBJDIR) -c -o $@ $<

$(OBJDIR)/mMolecule.o: $(SRCDIR)/mMolecule.f90 $(call obj,mStruc mGlobal mMath)
	$(FORT) $(FOPTS) -J$(OBJDIR) -c -o $@ $<

$(OBJDIR)/mCodeFlow.o: $(SRCDIR)/mCodeFlow.f90 $(OBJDIR)/mParser.o $(OBJDIR)/mStruc.o
	$(FORT) $(FOPTS) -J$(OBJDIR) -c -o $@ $<

$(OBJDIR)/mFix.o: $(SRCDIR)/mFix.f90 $(OBJDIR)/mGlobal.o
	$(FORT) $(FOPTS) -J$(OBJDIR) -c -o $@ $<

$(OBJDIR)/mMixingRule.o: $(SRCDIR)/mMixingRule.f90 $(call obj,mString mGlobal)
	$(FORT) $(FOPTS) -J$(OBJDIR) -c -o $@ $<

$(OBJDIR)/mParser.o: $(SRCDIR)/mParser.f90 $(OBJDIR)/mGlobal.o
	$(FORT) $(FOPTS) -J$(OBJDIR) -c -o $@ $<

$(OBJDIR)/mPackmol.o: $(SRCDIR)/mPackmol.f90 $(call obj,mMolecule mStruc mBox)
	$(FORT) $(FOPTS) -J$(OBJDIR) -c -o $@ $<

$(OBJDIR)/mBox.o: $(SRCDIR)/mBox.f90 $(OBJDIR)/mGlobal.o
	$(FORT) $(FOPTS) -J$(OBJDIR) -c -o $@ $<

$(OBJDIR)/mStruc.o: $(SRCDIR)/mStruc.f90 $(OBJDIR)/mString.o
	$(FORT) $(FOPTS) -J$(OBJDIR) -c -o $@ $<

$(OBJDIR)/mString.o: $(SRCDIR)/mString.f90 $(OBJDIR)/mGlobal.o
	$(FORT) $(FOPTS) -J$(OBJDIR) -c -o $@ $<

$(OBJDIR)/mMath.o: $(SRCDIR)/mMath.f90 $(OBJDIR)/mGlobal.o
	$(FORT) $(FOPTS) -J$(OBJDIR) -c -o $@ $<

$(OBJDIR)/mGlobal.o: $(SRCDIR)/mGlobal.f90
	mkdir -p $(OBJDIR)
	$(FORT) $(FOPTS) -J$(OBJDIR) -c -o $@ $<

