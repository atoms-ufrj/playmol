FORT  = gfortran
FOPTS = -O3 -march=native -ffast-math -funroll-loops -fstrict-aliasing -cpp -Wunused

SRCDIR  = ./src
OBJDIR  = $(SRCDIR)/obj
BINDIR  = ./bin
DOCDIR  = ./doc
GTKDIR  = ./highlight

PACKMOL = ./lib

exec = playmol

src = $(addprefix $(SRCDIR)/, $(addsuffix .f90, $(1)))
SRC  = $(call src, mBox mFix mMolecule mParser mString mCodeFlow mGlobal \
                   mPackmol mPlaymol mStruc playmol)
AUX  = $(call src, write_emdee write_lammps write_lammpstrj write_summary write_internals)
OBJ  = $(patsubst $(SRCDIR)/%.f90,$(OBJDIR)/%.o,$(SRC))

all: $(BINDIR)/$(exec)

.PHONY: install clean doc

doc:
	cd doc && doxygen

clean:
	rm -rf $(OBJDIR)
	rm -rf $(BINDIR)
	rm -rf $(DOCDIR)/html/
	cd $(PACKMOL) && make clean

install:
	cp -f $(BINDIR)/$(exec) /usr/local/bin
	cp -f $(PACKMOL)/packmol/packmol /usr/local/bin
	sh $(GTKDIR)/install.sh

$(BINDIR)/$(exec): $(OBJ) $(PACKMOL)/libpackmol.a
	mkdir -p $(BINDIR)
	$(FORT) $(FOPTS) -J$(OBJDIR) -o $@ $^

$(PACKMOL)/libpackmol.a:
	cd $(PACKMOL) && make

$(OBJDIR)/playmol.o: $(SRCDIR)/playmol.f90 $(OBJDIR)/mPlaymol.o
	$(FORT) $(FOPTS) -J$(OBJDIR) -c -o $@ $<

$(OBJDIR)/mPlaymol.o: $(SRCDIR)/mPlaymol.f90 $(AUX) $(OBJDIR)/mCodeFlow.o $(OBJDIR)/mPackmol.o \
                      $(OBJDIR)/mFix.o $(OBJDIR)/mBox.o $(OBJDIR)/mMolecule.o
	$(FORT) $(FOPTS) -J$(OBJDIR) -c -o $@ $<

$(OBJDIR)/mMolecule.o: $(SRCDIR)/mMolecule.f90 $(OBJDIR)/mStruc.o $(OBJDIR)/mGlobal.o
	$(FORT) $(FOPTS) -J$(OBJDIR) -c -o $@ $<

$(OBJDIR)/mCodeFlow.o: $(SRCDIR)/mCodeFlow.f90 $(OBJDIR)/mParser.o $(OBJDIR)/mStruc.o
	$(FORT) $(FOPTS) -J$(OBJDIR) -c -o $@ $<

$(OBJDIR)/mFix.o: $(SRCDIR)/mFix.f90 $(OBJDIR)/mGlobal.o
	$(FORT) $(FOPTS) -J$(OBJDIR) -c -o $@ $<

$(OBJDIR)/mParser.o: $(SRCDIR)/mParser.f90 $(OBJDIR)/mGlobal.o
	$(FORT) $(FOPTS) -J$(OBJDIR) -c -o $@ $<

$(OBJDIR)/mPackmol.o: $(SRCDIR)/mPackmol.f90 $(OBJDIR)/mMolecule.o \
                      $(OBJDIR)/mStruc.o $(OBJDIR)/mBox.o
	$(FORT) $(FOPTS) -J$(OBJDIR) -c -o $@ $<

$(OBJDIR)/mBox.o: $(SRCDIR)/mBox.f90 $(OBJDIR)/mGlobal.o
	$(FORT) $(FOPTS) -J$(OBJDIR) -c -o $@ $<

$(OBJDIR)/mStruc.o: $(SRCDIR)/mStruc.f90 $(OBJDIR)/mString.o
	$(FORT) $(FOPTS) -J$(OBJDIR) -c -o $@ $<

$(OBJDIR)/mString.o: $(SRCDIR)/mString.f90 $(OBJDIR)/mGlobal.o
	$(FORT) $(FOPTS) -J$(OBJDIR) -c -o $@ $<

$(OBJDIR)/mGlobal.o: $(SRCDIR)/mGlobal.f90
	mkdir -p $(OBJDIR)
	$(FORT) $(FOPTS) -J$(OBJDIR) -c -o $@ $<

