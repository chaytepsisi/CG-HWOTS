CC      ?= gcc
CFLAGS  ?= -O3 -march=native -Wall -Wextra -Wno-unused-parameter
LDFLAGS ?=
LIBS     = -lgmp -lm -ldl

SRC = src/sha256.c src/bqf.c src/cghwots.c src/main.c
OBJ = $(SRC:.c=.o)
BIN = cghwots

all: $(BIN)

$(BIN): $(OBJ)
	$(CC) $(CFLAGS) -o $@ $(OBJ) $(LDFLAGS) $(LIBS)

%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

clean:
	rm -f $(OBJ) $(BIN)

.PHONY: all clean
