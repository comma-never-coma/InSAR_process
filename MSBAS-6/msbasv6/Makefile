CC = g++
CXX = g++
CXXFLAGS = -fopenmp -Wall -I/usr/include/gdal # -Werror
LDFLAGS = -L/usr/lib
LDLIBS = -lgdal -llapack # -lf77blas -lcblas -latlas 

all: CArea.o  CImage.o  CInterferogram.o  CParam.o  CSet.o  main.o
	$(CC) $(CXXFLAGS) -o msbas CArea.o  CImage.o  CInterferogram.o  CParam.o  CSet.o  main.o $(LDFLAGS) $(LDLIBS)
	$(RM) *.o

clean: 
	$(RM) *.o msbas
