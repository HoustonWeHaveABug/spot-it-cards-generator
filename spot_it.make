SPOT_IT_C_FLAGS=-O2 -Wall -Wextra -Waggregate-return -Wcast-align -Wcast-qual -Wconversion -Wformat=2 -Winline -Wlong-long -Wmissing-prototypes -Wmissing-declarations -Wnested-externs -Wno-import -Wpointer-arith -Wredundant-decls -Wreturn-type -Wshadow -Wstrict-prototypes -Wswitch -Wwrite-strings

spot_it: spot_it.o
	gcc -o spot_it spot_it.o

spot_it.o: spot_it.c spot_it.make
	gcc -c ${SPOT_IT_C_FLAGS} -o spot_it.o spot_it.c

clean:
	rm -f spot_it spot_it.o
