objects = precmod.o dynarray.o stringmod.o xdr.o xdrwrap.o trajectoryio.o tranalmod.o tranal_base.o 

utilities = bileldens diel diffus dryrun h-bonds latdiff order rdf restime rgyr sdf torsion torsion2 trtcf 

FC=gfortran
 
STATIC = ../lib/libtranal_base.a

default: gfortran_set clean compile_all mklib

compile_all: $(utilities) 

gfortran_set: lxdr
	$(eval FFLAGS_OLD = -O3 -ffixed-form -g -fcheck=all -fbacktrace -cpp $(LXDR))
	$(eval FFLAGS = -O3 -ffree-form -g -fcheck=all -fbacktrace -cpp $(LXDR))
	$(eval TARGET= ../bin )

lxdr:	
# check if libxdrfile library is compiled and can be used
ifneq ($(wildcard ../lib/libxdrfile.a),)
	$(eval LXDR=-DXDR)
	$(eval LIBS= -L../lib -lxdrfile)
else
        $(info "!!!   File libxdrfile.a NOT found !!!")
        $(info "!!!   No XTC / TRR support           ")
endif

nolxdr:
	$(info "!!!----------------------------------------------------------------------!!!")
	$(info "Building without TRR / XTC support !")
	$(info "!!!----------------------------------------------------------------------!!!")
clean: 
	rm -f  *.o *.mod $(utilities)

bileldens: $(objects) bileldens.o 
	$(FC) $(FFLAGS) -o bileldens $^ $(LIBS)
	/bin/cp bileldens $(TARGET)

diel: $(objects) diel.o
	$(FC) $(FFLAGS) -o diel $^ $(LIBS)
	/bin/cp diel $(TARGET)

diffus: $(objects) diffus.o
	$(FC) $(FFLAGS) -o diffus $^ $(LIBS)
	/bin/cp diffus $(TARGET)

dryrun: $(objects) dryrun.o
	$(FC) $(FFLAGS) -o dryrun $^ $(LIBS)
	/bin/cp dryrun $(TARGET)

h-bonds: $(objects) h-bonds.o
	$(FC) $(FFLAGS) -o h-bonds $^ $(LIBS)
	/bin/cp h-bonds $(TARGET)

latdiff: $(objects) latdiff.o
	$(FC) $(FFLAGS) -o latdiff $^ $(LIBS)
	/bin/cp latdiff $(TARGET)

order: $(objects) order.o
	$(FC) $(FFLAGS) -o order $^ $(LIBS)
	/bin/cp order $(TARGET)

rdf: $(objects) rdf.o
	$(FC) $(FFLAGS) -o rdf $^ $(LIBS)
	/bin/cp rdf $(TARGET)

rgyr: $(objects) rgyr.o
	$(FC) $(FFLAGS) -o rgyr $^ $(LIBS)
	/bin/cp rgyr $(TARGET)

restime: $(objects) restime.o
	$(FC) $(FFLAGS) -o restime $^ $(LIBS)
	/bin/cp restime $(TARGET)

sdf: $(objects) sdf.o
	$(FC) $(FFLAGS) -o sdf $^ $(LIBS)
	/bin/cp sdf $(TARGET)

torsion: $(objects) torsion.o
	$(FC) $(FFLAGS) -o torsion $^ $(LIBS)
	/bin/cp torsion $(TARGET)

torsion2: $(objects) torsion2.o
	$(FC) $(FFLAGS) -o torsion2 $^ $(LIBS)
	/bin/cp torsion2 $(TARGET)

trtcf: $(objects) trtcf.o
	$(FC) $(FFLAGS) -o trtcf $^ $(LIBS)
	/bin/cp trtcf $(TARGET)

mklib: $(objects)
	@ar rcs $(STATIC) $(objects)

%.o: %.f90
	$(FC) $(FFLAGS) -c $< $(LIBS) -o $@
%.o: %.f
	$(FC) $(FFLAGS_OLD) -c $< $(LIBS) -o $@


