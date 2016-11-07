CC = gcc
CFLAGS = -I/usr/include/python2.7 -Wall -fPIC -O3


Hamil_wrap.c:
	swig -python Hamil.i

Hamil: Hamil_wrap.o
	ld -shared Hamil.o Hamil_wrap.o -o _Hamil.so 

%.o : %.c
	$(CC) -c Hamil.c $^ $(CFLAGS)


clean:
	rm -f ./*.o *~ ./Hamil.py ./Hamil_wrap.c ./*.pyc ./*.so

