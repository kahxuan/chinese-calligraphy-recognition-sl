ekmi.o: ekmi.c
	gcc -Wall \
		-I/usr/local/homebrew/cellar/libpng/1.6.37/include \
		-c ekmi.c -o ekmi.o && \
	gcc -Wall \
		-L/usr/local/homebrew/cellar/libpng/1.6.37/lib ekmi.o \
		-lpng -o ekmi.o

ekmi.so: ekmi.c
	gcc -shared -Wall  \
		-I/usr/local/homebrew/cellar/libpng/1.6.37/include \
		-L/usr/local/homebrew/cellar/libpng/1.6.37/lib \
		-lpng -o ekmi.so -fPIC ekmi.c 

ekmi.dylib: ekmi.c
	gcc -dynamiclib -Wall  \
		-I/usr/local/homebrew/cellar/libpng/1.6.37/include \
		-L/usr/local/homebrew/cellar/libpng/1.6.37/lib \
		-lpng -o ekmi.dylib -fPIC ekmi.c 

clean:
	rm *.o && rm *.so && rm *.dylib


