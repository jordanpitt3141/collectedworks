CC = gcc
CFLAGS = -I/usr/include/python2.7 -Wall -fPIC -O3


CoupAvec_wrap.c:
	swig -python CoupAvec.i

CoupAvec: CoupAvec_wrap.o
	ld -shared CoupAvec.o CoupAvec_wrap.o -o _CoupAvec.so 

%.o : %.c
	$(CC) -c CoupAvec.c $^ $(CFLAGS)


clean:
	rm -f ./*.o *~ ./CoupAvec.py ./CoupAvec_wrap.c ./*.pyc ./*.so

