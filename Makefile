     CC = gcc

CFLAGS = -pedantic -Wall -g -O2

LDFLAGS = -lm
     RM = rm -f

OBJS = stochmap.o mbmath_sub.o my_getopt.o

all : stochmap

stochmap : $(OBJS)
	$(CC) $(CFLAGS) -o stochmap $(OBJS) $(LDFLAGS)

%.o: %.c %.h
	$(CC) $(CFLAGS) -c $<

clean :
	$(RM) $(OBJS)
