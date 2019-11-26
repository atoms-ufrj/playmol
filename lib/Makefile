all: packmol
	make -C packmol/ lib
	mv packmol/libpackmol.a ./
	make -C packmol/

clean:
	rm -f libpackmol.a

clean-all:
	rm -rf packmol/
	rm -f libpackmol.a

packmol:
	git clone https://github.com/atoms-ufrj/packmol

