GFORTRAN=gfortran
FFLAGS=
PROG=hyper
VERSION=0.1
SOURCES=mult.f prep.f utility.f
OBJS=${SOURCES:.f=.o}

all: hyper

hyper: ${SOURCES}
	${GFORTRAN} ${FFLAGS} -o hyper ${SOURCES}

dist:
	@mkdir ${PROG}-${VERSION}
	@cp ${SOURCES} ${PROG}-${VERSION}
	@cp Makefile ${PROG}-${VERSION}
	@tar cvfz ${PROG}-${VERSION}.tar.gz ${PROG}-${VERSION} 
	@rm -rf ${PROG}-${VERSION}

clean:
	@rm -f ${PROG} ${OBJS} core

