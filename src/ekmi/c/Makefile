main: main.c
	gcc -Wall \
		-I/usr/local/homebrew/cellar/libpng/1.6.37/include \
		-c main.c -o main.o && \
	gcc -Wall \
		-L/usr/local/homebrew/cellar/libpng/1.6.37/lib main.o \
		-lpng -o main.o


ekmi.so: ekmi.c
	gcc -shared -Wall -o ../ekmi.so -fPIC ekmi.c 

ekmi.dylib: ekmi.c
	gcc -dynamiclib -Wall -o ../ekmi.dylib -fPIC ekmi.c 

clean:
	rm *.o