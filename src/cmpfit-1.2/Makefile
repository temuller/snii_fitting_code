CC=gcc
CFLAGS= -O3
OFILES=mpfit.o
LIBFILE=libmpfit.a

all: $(LIBFILE) 

clean:
	rm -f $(OFILES)  $(LIBFILE)

mpfit.o: mpfit.c mpfit.h
	$(CC) -c -o $@ $< $(CFLAGS)

$(LIBFILE): $(OFILES)
	$(AR) r $@ $(OFILES)
