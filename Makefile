INCLUDES = -Is lib

CFLAGS = 
LFLAGS = 

FLAGS = -ocamlopt "ocamlopt.opt -inline 100 -ffast-math"

LIBS = 

OCAMLBUILD = ocamlbuild
OCAMLBUILD_DOC = $(OCAMLBUILD) $(INCLUDES)
OCAMLBUILD_LIB = $(OCAMLBUILD_DOC) $(CFLAGS) $(LFLAGS) $(FLAGS)
OCAMLBUILD_EXE = $(OCAMLBUILD_LIB) $(LIBS)

.PHONY: all
all: doc lib

.PHONY: lib
lib: 
	$(OCAMLBUILD_LIB) nbody.cmxa nbody.cma

.PHONY: doc
doc:
	$(OCAMLBUILD_DOC) nbody.docdir/index.html

.PHONY: clean
clean:
	$(OCAMLBUILD) -clean