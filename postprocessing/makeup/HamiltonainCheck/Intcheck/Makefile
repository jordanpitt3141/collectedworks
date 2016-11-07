CC = gcc
CFLAGS = -I/usr/include/python2.7 -Wall -fPIC -O3


Int_wrap.c:
	swig -python Int.i

Int: Int_wrap.o
	ld -shared Int.o Int_wrap.o -o _Int.so 

%.o : %.c
	$(CC) -c Int.c $^ $(CFLAGS)


clean:
	rm -f ./*.o *~ ./Int.py ./Int_wrap.c ./*.pyc ./*.so

