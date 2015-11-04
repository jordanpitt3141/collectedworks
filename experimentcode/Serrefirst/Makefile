CC = gcc
CFLAGS = -I/usr/include/python2.7 -Wall -fPIC -O3


Serre1_wrap.c:
	swig -python Serre1.i

Serre1: Serre1_wrap.o
	ld -shared Serre1.o Serre1_wrap.o -o _Serre1.so 

%.o : %.c
	$(CC) -c Serre1.c $^ $(CFLAGS)


clean:
	rm -f ./*.o *~ ./Serre1.py ./Serre1_wrap.c ./*.pyc ./*.so

