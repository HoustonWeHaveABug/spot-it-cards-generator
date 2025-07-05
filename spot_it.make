SPOT_IT_C_FLAGS=-c -O2 -std=c89 -Wpedantic -Wall -Wextra -Waggregate-return -Wcast-align -Wcast-qual -Wconversion -Wformat=2 -Winline -Wlong-long -Wmissing-prototypes -Wmissing-declarations -Wnested-externs -Wpointer-arith -Wredundant-decls -Wshadow -Wstrict-prototypes -Wwrite-strings -Wswitch-default -Wswitch-enum -Wbad-function-cast -Wstrict-overflow=5 -Wundef -Wlogical-op -Wfloat-equal -Wold-style-definition

spot_it: spot_it.o
	gcc -o spot_it spot_it.o

spot_it.o: spot_it.c spot_it.make
	gcc ${SPOT_IT_C_FLAGS} -o spot_it.o spot_it.c

clean:
	rm -f spot_it spot_it.o
