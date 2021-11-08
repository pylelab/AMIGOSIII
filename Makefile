CC=g++
CFLAGS=-O3
LDFLAGS=#-static

NaTorsion: NaTorsion.cpp NaTorsion.hpp PDBParser.hpp pstream.h GeometryTools.hpp
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}
