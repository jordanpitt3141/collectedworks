CC = gcc
CFLAGS = -I/usr/include/python2.7 -Wall -fPIC -O3


Serre2_wrap.c:
	swig -python Serre2.i

Serre2: Serre2_wrap.o
	ld -shared Serre2.o Serre2_wrap.o -o _Serre2.so 

%.o : %.c
	$(CC) -c Serre2.c $^ $(CFLAGS)


clean:
	rm -f ./*.o *~ ./Serre2.py ./Serre2_wrap.c ./*.pyc ./*.so

