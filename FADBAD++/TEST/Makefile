CC = g++
CFLAGS = -I..
OBJ = TestFADBAD.o TestTAD.o TestAll.o

.cpp.o:
	$(CC) $(CFLAGS) -c $<

TestAll : $(OBJ)
	$(CC) -o $@ $(OBJ)

all: TestAll

clean:
	-rm *.o TestAll