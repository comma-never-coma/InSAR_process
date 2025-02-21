CC = g++
CXX = g++
CXXFLAGS = -fopenmp -Wall -I/usr/include/gdal # -Werror
LDFLAGS = -L/usr/lib
LDLIBS = -lgdal

all: msbas_extract.o
	$(CC) $(CXXFLAGS) -o msbas_extract msbas_extract.o $(LDFLAGS) $(LDLIBS)
	$(RM) *.o

clean: 
	$(RM) *.o msbas_extract
