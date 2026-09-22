# Generic makefile for projects with SRC, BIN, OBJ 
# hierarchy

# setting the shell environment
SHELL = /bin/bash

# setting directories
SRCDIR = SRC
INCDIR = INC
OBJDIR = OBJ
BINDIR = BIN

# compilation environment
CC = gcc
CFLAGS = -O3 -std=c99
#CFLAGS = -g -O0
LFLAGS = -lm -llapack $(shell pkg-config --libs fftw3) -lblas
INCDIRS = -I${INCDIR} $(shell pkg-config --cflags fftw3)


SOURCES  = $(wildcard $(SRCDIR)/*.c)
INCS  = $(wildcard $(INCDIR)/*.h)
OBJECTS  = $(SOURCES:$(SRCDIR)/%.c=$(OBJDIR)/%.o)
TARGET = $(BINDIR)/trilat.x 
rm       = rm -f

$(TARGET): $(OBJECTS)
	$(CC) $(OBJECTS) $(CFLAGS) $(LFLAGS) $(INCDIRS) -o $(TARGET)

$(OBJDIR)/%.o : $(SRCDIR)/%.c $(INCS)
	$(CC) -c $(CFLAGS) $(LFLAGS) $(INCDIRS) $< -o $@

clean:
	$(rm) $(OBJDIR)/*.o $(BINDIR)/*.x
