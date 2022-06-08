CC=mpicc
CFLAGS=-Wall -Wextra -pedantic -O3 -g -pthread
#CFLAGS= -pedantic -O3 -g -pthread
LIBS=-lm -lssl -lcrypto


all: attack rainbow

###attack
attack: attack.o
	mkdir output
	gcc -o $@ $^ $(LIBS) 

### generation
rainbow: table.o rainbow.o crypto.o io.o
	mkdir datas
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS) $(LDFLAGS)

table.o: rainbow_scale.h table.c table.h
	$(CC) $(CFLAGS) -o table.o -c table.c

crypto.o: crypto.c crypto.h
	$(CC) $(CFLAGS) -o crypto.o -c crypto.c

rainbow.o: rainbow_scale.c rainbow_scale.h
	$(CC) $(CFLAGS) -o rainbow.o -c rainbow_scale.c

io.o: io.c io.h rainbow_scale.h
	$(CC) $(CFLAGS) -o io.o -c io.c

attack.o: attack.c
	gcc -c attack.c 

clean:
	rm -rf *.o 
	rm -rf datas
	rm -rf output
	rm -rf rainbow
	rm -rf attack
