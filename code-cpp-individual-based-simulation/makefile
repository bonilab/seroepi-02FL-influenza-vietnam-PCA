# driver: driver.o rkf.o
# 	lcc  driver.o rkf.o -o driver -lm
# makefile
# Tabs *must* be used for the indentations below;
# spaces cause make syntax errors.

CC=g++
#CFLAGS=-fast -xO4 -xdepend -xarch=v8plusa -xprefetch -xvector -xunroll=8 -fsimple=2 -xsafe=mem
LIBS=-lm
GSLLIBS=-lgsl -lgslcblas 

runsim:		
		$(CC) -O3 -c -o person.o person.cpp $(LIBS) $(GSLLIBS)
		$(CC) -O3 -c -o readdata.o readdata.cpp $(LIBS) $(GSLLIBS)
		$(CC) -O3 -c -o main_sim.o main_sim.cpp $(LIBS) $(GSLLIBS)
		$(CC) $(CFLAGS) -o runsim main_sim.o readdata.o person.o $(LIBS) $(GSLLIBS) 
		

clean:
		rm -f *.o core runsim *~




