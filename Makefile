FORT  = gfortran
FOPTS = -O3 -march=native -ffast-math -funroll-loops -fstrict-aliasing -cpp -Wunused

SRCDIR  = ./src
OBJDIR  = $(SRCDIR)/obj
BINDIR  = ./bin
DOCDIR  = ./doc

PACKMOL = ./lib

exec = playmol
src  = ./src

all: $(BINDIR)/$(exec)

.PHONY: install doc clean

doc:
	doxygen $(DOCDIR)/Doxyfile

clean:
	rm -rf $(OBJDIR)
	rm -rf $(BINDIR)
	rm -rf $(DOCDIR)/html/
	rm -rf $(DOCDIR)/latex/

install:
	cp -f $(BINDIR)/$(exec) /usr/local/bin

$(BINDIR)/$(exec): $(OBJDIR)/playmol.o $(OBJDIR)/mData.o $(OBJDIR)/mStruc.o $(OBJDIR)/mString.o $(OBJDIR)/mGlobal.o
	mkdir -p $(BINDIR)
	$(FORT) $(FOPTS) -J$(OBJDIR) -L$(PACKMOL) -lpackmol -o $@ $^

$(OBJDIR)/playmol.o: $(SRCDIR)/playmol.f90 $(OBJDIR)/mData.o $(OBJDIR)/mStruc.o $(OBJDIR)/mString.o $(OBJDIR)/mGlobal.o
	$(FORT) $(FOPTS) -J$(OBJDIR) -c -o $@ $<

$(OBJDIR)/mData.o: $(SRCDIR)/mData.f90 $(OBJDIR)/mStruc.o
	$(FORT) $(FOPTS) -J$(OBJDIR) -c -o $@ $<

$(OBJDIR)/mStruc.o: $(SRCDIR)/mStruc.f90 $(OBJDIR)/mString.o $(OBJDIR)/mGlobal.o
	$(FORT) $(FOPTS) -J$(OBJDIR) -c -o $@ $<

$(OBJDIR)/mString.o: $(SRCDIR)/mString.f90 $(OBJDIR)/mGlobal.o
	$(FORT) $(FOPTS) -J$(OBJDIR) -c -o $@ $<

$(OBJDIR)/mGlobal.o: $(SRCDIR)/mGlobal.f90
	mkdir -p $(OBJDIR)
	$(FORT) $(FOPTS) -J$(OBJDIR) -c -o $@ $<

