main: main.c
	gcc -Wall \
		-I/opt/homebrew/cellar/libpng/1.6.37/include \
		-c main.c -o out.o && \
	gcc -Wall \
		-L/opt/homebrew/cellar/libpng/1.6.37/lib out.o \
		-lgsl -lpng -o out.e

clean: *.o *.e
	rm *.o *.e