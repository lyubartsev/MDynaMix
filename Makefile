default:
	make -f Makefile.gfortran

linux:	
	make -f Makefile.gfortran
f77:	
	make -f Makefile.f77
pgi:	
	make -f Makefile.pgi  
intel:	
	make -f Makefile.ifort
static:
	make -f Makefile.static
mpi:
	make -f Makefile.mpi
clean:
	rm -f  *.o
