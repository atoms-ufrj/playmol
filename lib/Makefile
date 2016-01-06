VERSION = 15.287

all:
	tar -zxvf packmol-$(VERSION).tar.gz
	cd packmol && sed -i 's/\r//g' *.f90 && patch < ../packmol.patch && make && make lib
	mv packmol/libpackmol.a ./

clean:
	rm -rf packmol/
	rm -f libpackmol.a

