CFLAGS=-O2 -Wall -Wextra -g
LDFLAGS=-lm

TARGET=zhunt
SOURCES=zhunt3.c antisyn.c delta_linking.c

all: $(TARGET)

$(TARGET): $(SOURCES)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

clean:
	rm -f $(TARGET)

.PHONY: clean
