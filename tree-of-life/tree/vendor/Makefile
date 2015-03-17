CC=gcc
CFLAGS=-I.
DEPS = rmq.h
OBJ = rmqtest.o

%.o: %.c $(DEPS)
		$(CC) -c -o $@ $< $(CFLAGS)

rmqtest: $(OBJ)
		gcc -o $@ $^ $(CFLAGS)

