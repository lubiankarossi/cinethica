CC=gcc
#CFLAGS=-W -Wall -ansi -pedantic -lm -std=c99
CFLAGS= -lm -std=c99
EXEC=cinethica
OBJ=cinetica.o model.o helper.o inp.o
all: $(EXEC)
cinethica: cinetica.o model.o helper.o inp.o
	$(CC) -o $@ $^ $(CFLAGS)
%.o: %.c
	$(CC) -o $@ -c $< $(CFLAGS)
main.o: reatimetro.c model.o pile.o
	$(CC) -o $@ -c $< $(CFLAGS)
.PHONY: clean mrproper

clean: 
	rm -rf *.o *~ $(EXEC)

mrproper: clean
	rm -rf $(EXEC)

