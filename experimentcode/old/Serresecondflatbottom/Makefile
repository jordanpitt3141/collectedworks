CC = gcc
CFLAGS = -I/usr/include/python2.7 -Wall -fPIC -O3


Serre2dc_wrap.c:
	swig -python Serre2dc.i

Serre2dc: Serre2dc_wrap.o
	ld -shared Serre2dc.o Serre2dc_wrap.o -o _Serre2dc.so 

%.o : %.c
	$(CC) -c Serre2dc.c $^ $(CFLAGS)


clean:
	rm -f ./*.o *~ ./Serre2dc.py ./Serre2dc_wrap.c ./*.pyc ./*.so

