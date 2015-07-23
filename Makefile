CC=gcc
CFLAGS=-Wall -Wextra
LDFLAGS=-lrt
SOURCES=bam-errorrate.c
OBJECTS=$(SOURCES:.c=.o)
EXECUTABLE=bam-errorrate

ifeq ($(DEBUG), 1)
	CFLAGS+=-g
else
	CFLAGS+=-O3
endif

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.c:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -vf $(EXECUTABLE) $(OBJECTS)
