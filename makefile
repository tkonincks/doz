EXEC=doz.x
GC=gfortran-4.7
GOPT=-O3 -Wall

$(EXEC): main-doz.o read.o struct.o eigen.o fft.o fileman.o func.o 
	$(GC) $(GOPT) -o $(EXEC) main-doz.o read.o struct.o eigen.o fft.o fileman.o func.o

main-doz.o: main-doz.f95
	$(GC) $(GOPT) -c main-doz.f95

read.o: read.f95
	$(GC) $(GOPT) -c read.f95

struct.o: struct.f95
	$(GC) $(GOPT) -c struct.f95

eigen.o: eigen.f95
	$(GC) $(GOPT) -c eigen.f95

fft.o: fft.f
	$(GC) $(GOPT) -c fft.f

fileman.o: fileman.f95
	$(GC) $(GOPT) -c fileman.f95

func.o: func.f95
	$(GC) $(GOPT) -c func.f95

clean:
	rm -rf *.o *.dat

mrproper: clean
	rm -rf doz.x *.o *.dat
