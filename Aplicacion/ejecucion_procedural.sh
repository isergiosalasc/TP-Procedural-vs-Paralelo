#!/bin/bash
# Variables
muestras=5
resta=0
start_time=$(date +%s)

for i in $( seq 1 $muestras); do
	./TP2_procedural.out
done

finish_time=$(date +%s)
resta=$((finish_time - start_time))
tiempo_total=$( calc $resta/$muestras)

echo
echo ">	El tiempo de demora total para las $muestras muestras es: $resta segundos."
echo ">	El timepo promedio por cada muestra es:" $tiempo_total "segundos."