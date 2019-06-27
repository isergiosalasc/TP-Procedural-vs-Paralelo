#!/bin/bash
# Variables
muestras=10
hilos=4

for i in $( seq 1 $muestras); do
	./TP2_paralelo_2.out $hilos
done