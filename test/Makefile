PROG := tsignal
SHELL := /bin/sh
CC := gcc
CXX := g++

COMPILER_OPTIONS := -m64 -Wall -Wextra -Wshadow -Werror -pedantic -Iinclude
CFLAGS := -std=c99 $(COMPILER_OPTIONS)
CXXFLAGS := -std=c++11 -Weffc++ $(COMPILER_OPTIONS)
LDFLAGS := -Wl,--no-as-needed -lm

DEBUGFLAGS := -g -O0 -D _DEBUG
RELEASEFLAGS := -O2 -D NDEBUG

ifeq ($(PREFIX),)
	PREFIX := /usr/local
endif

TARGET  := test-$(PROG).x
SOURCES := test_$(PROG).cpp
LIBS 	:= -l$(PROG) -lfftw3 -lunittest
HEADERS := $(PREFIX)/include/$(PROG).h
COMMON	:=
OBJECTS := $(SOURCES:.cpp=.o)

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) $(FLAGS) $(CXXFLAGS) $(DEBUGFLAGS) -o $(TARGET) $(OBJECTS) $(LIBS)

clean:
	-rm -f $(TARGET)
	-rm -f $(OBJECTS)

%.o: %.c $(HEADERS) $(COMMON)
	$(CXX) $(CXXFLAGS) $(DEBUGFLAGS) -c -o $@ $<

.PHONY: all clean
