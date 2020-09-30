objects = aver.o cpu_intel.o input.o forces.o main.o mdee.o mdstep.o \
       restart.o scalar.o service.o setup.o util.o getcpu.o

   FC = f77
   FFLAGS = 
   LDFLAGS = 
   TARGET = md

default: $(objects) 
	$(FC) $(LDFLAGS) -o $(TARGET) $(objects)
   $(objects) : prcm.h dimpar.h

clean: 
	rm -f $(TARGET) *.o
