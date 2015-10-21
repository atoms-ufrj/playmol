VERSION = 14.072

all:
	tar -zxvf packmol-$(VERSION).tar.gz
	cd packmol && patch < ../packmol.patch && make && make lib
	mv packmol/libpackmol.a ./

clean:
	rm -rf packmol/
	rm -f libpackmol.a

