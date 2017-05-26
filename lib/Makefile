all: packmol
	make -C packmol/ lib
	mv packmol/libpackmol.a ./

clean:
	rm -f libpackmol.a

cleanall:
	rm -rf packmol/
	rm -f libpackmol.a

packmol:
	git clone https://github.com/atoms-ufrj/packmol

