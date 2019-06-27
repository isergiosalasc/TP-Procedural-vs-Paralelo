# Makefile del TP2

#VARIABLES

#REGLAS
make:
	gcc -fopenmp -g -pg -Wall -Werror -pedantic TP2_interfaz.c -o TP2_interfaz.out
	gcc -fopenmp -g -pg -Wall -Werror -pedantic TP2_procedural.c -o TP2_procedural.out -lm -std=gnu99
	gcc -fopenmp -g -pg -Wall -Werror -pedantic TP2_paralelo.c -o TP2_paralelo.out -lm -std=gnu99
	gcc -fopenmp -g -pg -Wall -Werror -pedantic TP2_paralelo_2.c -o TP2_paralelo_2.out -lm -std=gnu99
clean:
	rm -rf *.out
	rm -rf *~
	rm -rf resultados_adpc
ejecutar:
	./TP2_procedural