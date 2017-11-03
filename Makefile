CC = gcc
CFLAGS = -O3 -Wall -Werror -pedantic -fopenmp
LDFLAGS = -lm
TARGET = remd
OBJECTS = config.o init.o integrator.o potential.o main.o

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) -o $@ $(LDFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) -c $<

clean:
	rm -f $(OBJECTS)
