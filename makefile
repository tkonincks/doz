EXEC=doz.x
#GC=gfortran-4.7
GC=ifort

ifeq ($(GC),ifort)

GOPT=-mtune=generic -O2 -xHost -ipo

# -msse4a

$(EXEC): main-doz.o read.o struct.o eigen.o dyn.o fft.o fileman.o func.o 
	$(GC) -o $(EXEC) main-doz.o read.o struct.o eigen.o dyn.o fft.o fileman.o func.o

main-doz.o: main-doz.f90
	$(GC) -c $(GOPT) main-doz.f90
                
read.o: read.f90
	$(GC) -c $(GOPT) read.f90
                
struct.o: struct.f90
	$(GC) -c $(GOPT) struct.f90
                
eigen.o: eigen.f90
	$(GC) -c $(GOPT) eigen.f90
                
dyn.o: dyn.f90
	$(GC) -c $(GOPT) dyn.f90

fft.o: fft.f    
	$(GC) -c $(GOPT) fft.f
                
fileman.o: fileman.f90
	$(GC) -c $(GOPT) fileman.f90
                
func.o: func.f90
	$(GC) -c $(GOPT) func.f90

else

GOPT=-Wall -O2

$(EXEC): main-doz.o read.o struct.o eigen.o dyn.o fft.o fileman.o func.o 
	$(GC) -o $(EXEC) main-doz.o read.o struct.o eigen.o dyn.o fft.o fileman.o func.o

main-doz.o: main-doz.f90
	$(GC) $(GOPT) -c main-doz.f90

read.o: read.f90
	$(GC) $(GOPT) -c read.f90

struct.o: struct.f90
	$(GC) $(GOPT) -c struct.f90

eigen.o: eigen.f90
	$(GC) $(GOPT) -c eigen.f90
                
dyn.o: dyn.f90
	$(GC) -c $(GOPT) dyn.f90

fft.o: fft.f
	$(GC) $(GOPT) -c fft.f

fileman.o: fileman.f90
	$(GC) $(GOPT) -c fileman.f90

func.o: func.f90
	$(GC) $(GOPT) -c func.f90

endif

clean:
	rm -rf *.o *.dat

mrproper: clean
	rm -rf doz.x *.o *.dat
