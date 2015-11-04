CC = gcc
CFLAGS = -I/usr/include/python2.7 -Wall -fPIC -O3


Serre2GR_wrap.c:
	swig -python Serre2GR.i

Serre2GR: Serre2GR_wrap.o
	ld -shared Serre2GR.o Serre2GR_wrap.o -o _Serre2GR.so 

%.o : %.c
	$(CC) -c Serre2GR.c $^ $(CFLAGS)


clean:
	rm -f ./*.o *~ ./Serre2GR.py ./Serre2GR_wrap.c ./*.pyc ./*.so

