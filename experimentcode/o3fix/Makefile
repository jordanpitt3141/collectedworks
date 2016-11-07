CC = gcc
CFLAGS = -I/usr/include/python2.7 -Wall -fPIC -O3


Serre3_wrap.c:
	swig -python Serre3.i

Serre3: Serre3_wrap.o
	ld -shared Serre3.o Serre3_wrap.o -o _Serre3.so 

%.o : %.c
	$(CC) -c Serre3.c $^ $(CFLAGS)


clean:
	rm -f ./*.o *~ ./Serre3.py ./Serre3_wrap.c ./*.pyc ./*.so

