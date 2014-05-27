FORT  = gfortran
FOPTS = -O3 -march=native -ffast-math -funroll-loops -fstrict-aliasing -cpp -Wunused

SRCDIR  = ./src
OBJDIR  = $(SRCDIR)/obj
DOCDIR  = ./doc

exec = ./bin/playmol
src  = ./src

all: $(exec)

.PHONY: install doc clean

doc:
	doxygen $(DOCDIR)/Doxyfile

clean:
	rm -rf $(OBJDIR)
	rm -rf $(exec)
	rm -rf $(DOCDIR)/html/
	rm -rf $(DOCDIR)/latex/

install:
	cp -f $(exec) /usr/local/bin

$(exec): $(SRCDIR)/playmol.f90 $(OBJDIR)/mData.o $(OBJDIR)/mStruc.o $(OBJDIR)/mString.o $(OBJDIR)/mGlobal.o
	$(FORT) $(FOPTS) -J$(OBJDIR) -o $@ $^

$(OBJDIR)/mData.o: $(SRCDIR)/mData.f90 $(OBJDIR)/mStruc.o
	$(FORT) $(FOPTS) -J$(OBJDIR) -c -o $@ $<

$(OBJDIR)/mStruc.o: $(SRCDIR)/mStruc.f90 $(OBJDIR)/mString.o $(OBJDIR)/mGlobal.o
	$(FORT) $(FOPTS) -J$(OBJDIR) -c -o $@ $<

$(OBJDIR)/mString.o: $(SRCDIR)/mString.f90 $(OBJDIR)/mGlobal.o
	$(FORT) $(FOPTS) -J$(OBJDIR) -c -o $@ $<

$(OBJDIR)/mGlobal.o: $(SRCDIR)/mGlobal.f90
	mkdir -p $(OBJDIR)
	$(FORT) $(FOPTS) -J$(OBJDIR) -c -o $@ $<

