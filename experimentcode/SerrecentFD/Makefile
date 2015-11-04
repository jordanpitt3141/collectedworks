CC = gcc
CFLAGS = -I/usr/include/python2.7 -Wall -fPIC -O3


Serre2FDC_wrap.c:
	swig -python Serre2FDC.i

Serre2FDC: Serre2FDC_wrap.o
	ld -shared Serre2FDC.o Serre2FDC_wrap.o -o _Serre2FDC.so 

%.o : %.c
	$(CC) -c Serre2FDC.c $^ $(CFLAGS)


clean:
	rm -f ./*.o *~ ./Serre2FDC.py ./Serre2FDC_wrap.c ./*.pyc ./*.so

