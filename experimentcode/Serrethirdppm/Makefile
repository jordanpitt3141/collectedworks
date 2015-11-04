CC = gcc
CFLAGS = -I/usr/include/python2.7 -Wall -fPIC -O3


Serre3ppm_wrap.c:
	swig -python Serre3ppm.i

Serre3ppm: Serre3ppm_wrap.o
	ld -shared Serre3ppm.o Serre3ppm_wrap.o -o _Serre3ppm.so 

%.o : %.c
	$(CC) -c Serre3ppm.c $^ $(CFLAGS)


clean:
	rm -f ./*.o *~ ./Serre3ppm.py ./Serre3ppm_wrap.c ./*.pyc ./*.so

