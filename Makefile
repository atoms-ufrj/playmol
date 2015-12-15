FORT  = gfortran
FOPTS = -O3 -march=native -ffast-math -funroll-loops -fstrict-aliasing -cpp -Wunused

SRCDIR  = ./src
OBJDIR  = $(SRCDIR)/obj
BINDIR  = ./bin
DOCDIR  = ./doc
GTKDIR  = ./highlight

PACKMOL = ./lib

exec = playmol
src  = ./src

all: $(BINDIR)/$(exec)

.PHONY: install clean doc

doc:
	doxygen $(DOCDIR)/Doxyfile

clean:
	rm -rf $(OBJDIR)
	rm -rf $(BINDIR)
	rm -rf $(DOCDIR)/html/
	rm -rf $(DOCDIR)/latex/
	cd $(PACKMOL) && make clean

install:
	cp -f $(BINDIR)/$(exec) /usr/local/bin
	cp -f $(PACKMOL)/packmol/packmol /usr/local/bin
	sh $(GTKDIR)/install.sh

$(BINDIR)/$(exec): $(OBJDIR)/playmol.o $(OBJDIR)/mPlaymol.o $(OBJDIR)/mStruc.o  \
                   $(OBJDIR)/mPackmol.o $(OBJDIR)/mBox.o $(OBJDIR)/mString.o    \
                   $(OBJDIR)/mAlign.o $(OBJDIR)/mParser.o $(OBJDIR)/mGlobal.o   \
                   $(PACKMOL)/libpackmol.a
	mkdir -p $(BINDIR)
	$(FORT) $(FOPTS) -J$(OBJDIR) -o $@ $^

$(PACKMOL)/libpackmol.a:
	cd $(PACKMOL) && make

$(OBJDIR)/playmol.o: $(SRCDIR)/playmol.f90 $(OBJDIR)/mPlaymol.o
	$(FORT) $(FOPTS) -J$(OBJDIR) -c -o $@ $<

$(OBJDIR)/mPlaymol.o: $(SRCDIR)/mPlaymol.f90 $(OBJDIR)/mPackmol.o $(OBJDIR)/mBox.o $(OBJDIR)/mAlign.o $(OBJDIR)/mParser.o
	$(FORT) $(FOPTS) -J$(OBJDIR) -c -o $@ $<

$(OBJDIR)/mPackmol.o: $(SRCDIR)/mPackmol.f90 $(OBJDIR)/mStruc.o $(OBJDIR)/mBox.o
	$(FORT) $(FOPTS) -J$(OBJDIR) -c -o $@ $<

$(OBJDIR)/mStruc.o: $(SRCDIR)/mStruc.f90 $(OBJDIR)/mString.o
	$(FORT) $(FOPTS) -J$(OBJDIR) -c -o $@ $<

$(OBJDIR)/mBox.o: $(SRCDIR)/mBox.f90 $(OBJDIR)/mGlobal.o
	$(FORT) $(FOPTS) -J$(OBJDIR) -c -o $@ $<

$(OBJDIR)/mAlign.o: $(SRCDIR)/mAlign.f90 $(OBJDIR)/mGlobal.o
	$(FORT) $(FOPTS) -J$(OBJDIR) -c -o $@ $<

$(OBJDIR)/mParser.o: $(SRCDIR)/mParser.f90 $(OBJDIR)/mString.o
	$(FORT) $(FOPTS) -J$(OBJDIR) -c -o $@ $<

$(OBJDIR)/mString.o: $(SRCDIR)/mString.f90 $(OBJDIR)/mGlobal.o
	$(FORT) $(FOPTS) -J$(OBJDIR) -c -o $@ $<

$(OBJDIR)/mGlobal.o: $(SRCDIR)/mGlobal.f90
	mkdir -p $(OBJDIR)
	$(FORT) $(FOPTS) -J$(OBJDIR) -c -o $@ $<

