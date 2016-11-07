CC = gcc
CFLAGS = -I/usr/include/python2.7 -Wall -fPIC -O3


Serre3fem_wrap.c:
	swig -python Serre3fem.i

Serre3fem: Serre3fem_wrap.o
	ld -shared Serre3fem.o Serre3fem_wrap.o -o _Serre3fem.so 

%.o : %.c
	$(CC) -c Serre3fem.c $^ $(CFLAGS)


clean:
	rm -f ./*.o *~ ./Serre3fem.py ./Serre3fem_wrap.c ./*.pyc ./*.so

