PROG := signal
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

SRCDIR	:= src
INCDIR	:= include
OBJDIR 	:= obj
BINDIR	:= bin
LIBDIR	:= lib

HEADER := $(PROG).h
SOURCES := $(shell echo $(SRCDIR)/*.cpp)
HEADERS := $(shell echo $(INCDIR)/*.h)
COMMON  :=
OBJECTS := $(patsubst $(SRCDIR)/%,$(OBJDIR)/%,$(SOURCES:.cpp=.o))

TARGET  := $(LIBDIR)/lib$(PROG).a

all: $(TARGET)

$(TARGET): $(OBJECTS)
	ar ru $@ $^
	ranlib $@

test:
	cd test && $(MAKE)

release: $(SOURCES) $(HEADERS) $(COMMON)
	$(CXX) $(FLAGS) $(CXXFLAGS) $(RELEASEFLAGS) -o $(TARGET) $(SOURCES)

install: $(TARGET)
	install -d $(PREFIX)/lib/
	install -m 644 $(TARGET) $(PREFIX)/lib/
	install -d $(PREFIX)/include/
	install -m 644 include/$(HEADER) $(PREFIX)/include/

zip:
	-zip $(PROG).zip $(HEADERS) $(SOURCES) Makefile

clean:
	-rm -f $(TARGET) $(OBJECTS) $(PROG).zip
	-rm -f $(PREFIX)/lib/$(TARGET) $(PREFIX)/include/$(HEADER)

remake: clean all

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp $(HEADERS) $(COMMON)
	$(CXX) $(CXXFLAGS) $(DEBUGFLAGS) -c -o $@ $<

.PHONY : all release install zip clean test remake
