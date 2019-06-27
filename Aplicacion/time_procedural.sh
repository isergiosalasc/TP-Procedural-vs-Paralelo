#!/bin/bash
# Variables
muestras=5
hilos=4

for i in $( seq 1 $muestras); do
	./TP2_procedural.out
done