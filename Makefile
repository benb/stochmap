     CC = gcc

CFLAGS = -pedantic -Wall -g -O2

LDFLAGS = -lm
     RM = rm -f

OBJS = stochmap.o mbmath_sub.o my_getopt.o
OBJSL = libstochmap.o mbmath_sub.o my_getopt.o

all : stochmap

stochmap : $(OBJS)
	$(CC) $(CFLAGS) -o stochmap $(OBJS) $(LDFLAGS)

libstochmap.o: stochmap.c stochmap.h
	$(CC) $(CFLAGS) -fPIC -c -DBUILDLIBRARY -o $@ $<


%.o: %.c %.h
	$(CC) $(CFLAGS) -fPIC -c $<

lib: libstochmap.so.1.0.1
libstochmap.so.1.0.1: $(OBJSL)
	gcc -shared -Wl,-soname,libstochmap.so.1 -o libstochmap.so.1.0.1 $(OBJSL) 
clean :
	$(RM) $(OBJS)
