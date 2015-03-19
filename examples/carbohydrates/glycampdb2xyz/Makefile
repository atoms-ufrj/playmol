FORT  = gfortran
FOPTS = -O3 -march=native -ffast-math -funroll-loops -fstrict-aliasing -cpp -Wunused

EXEC = glycampdb2xyz

all: $(EXEC)

.PHONY: install clean

clean:
	rm -rf *.o
	rm -rf *.mod
	rm -rf $(EXEC)

install:
	cp -f $(EXEC) /usr/local/bin

$(EXEC): $(EXEC).o mString.o mGlobal.o
	$(FORT) $(FOPTS) -o $@ $^

$(EXEC).o: $(EXEC).f90 mString.o mGlobal.o
	$(FORT) $(FOPTS) -c -o $@ $<

mString.o: mString.f90 mGlobal.o
	$(FORT) $(FOPTS) -c -o $@ $<

mGlobal.o: mGlobal.f90
	$(FORT) $(FOPTS) -c -o $@ $<

