#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <unistd.h>
#include <omp.h>
#include <string.h>

#define GATES 500
#define PULSOS 100*8 /* 100 pulsos por 8 grado de acimut considerados */
#define name_size 100

/* Declaracion de funciones utilizadas en aplicacion */
void gate_por_canal(char *f_name, float canal_v[GATES][PULSOS], float canal_h[GATES][PULSOS]);
void autocorrelacion(float canal_v[GATES][PULSOS], float canal_h[GATES][PULSOS], float p_autoc_v[], float p_autoc_h[]);
int abrir_binario_resultados(char *f_name);
void ver_resultados(float autoc_v[], float autoc_h[]);
int convertir_cadena(char *pc);
void ver_estadisticas(double tp_est, int im);


/* Funcion principal */
int main(int argc, char const *argv[]) {

	/* Nombres de ficheros utilizados*/
	char f_pulsos[name_size]="pulsos.iq";
	char fout_name[name_size]="resultados_adpc"; /* resultados de autocorrelacion discreta por pulsos de cada canal*/

	
	/*Variables para contar tiempo de ejecución*/
	double total_i, total_f;

	/* Matrices de variables para las estructura gate-pulso por canal*/
	float canal_v[GATES][PULSOS];
	float canal_h[GATES][PULSOS];
	
	/* Vectores con los resultados de la autocorrelacion discreta por pulsos de cada canal*/
	float autoc_v[GATES];
	float autoc_h[GATES];

	printf("\nEl numero de pulsos que tiene el archivo %s es: %d \n",f_pulsos, PULSOS);
	printf("\nLas componetes de rango asociadas a cada gate son: %d \n", GATES);
	printf("\nSe utlizan estructuras de datos del tipo matriz[GATES][PULSO] = canal[%d][%d] \n", GATES, PULSOS);

	int opcion =0;
	int i,k=0;
	float t_est=0; /* tiempo estadistico segun numero de iteraciones por algoritmo*/
	int im=0; /*iterador de muestras*/
	int iteraciones=1;
	char cadena[256];
	/* Variables de archivos.*/
	FILE *f_out;
	int cant_gates= (int)GATES;
	while(opcion!=6)
	{
		printf( "\n   >>>_			MENU 			_<<<\n");
        printf( "\n   1. Cambiar numero de iteraciones por algoritmo.");
        printf( "\n   2. Resolver Algoritmo de autocerrelacion.");
        printf( "\n   3. Ver estadisticas de algoritmos.");
        printf( "\n   4. Mostrar resultados desde archivo binario %s.", fout_name);
        printf( "\n   5. Mostrar resultados desde estructura autocorrelacion.");
        printf( "\n   6. Salir." );
        printf( "\n\n   Introduzca opcion (1-6): ");

        scanf( "%d", &opcion );
		switch(opcion)
		{
			case 1:	memset(cadena, '\0', 256);
	                printf( "\n   > Introduzca el numero de iteraciones deseado: ");
	                scanf( "%s", cadena);
	                iteraciones = convertir_cadena(cadena);
	                break;

			case 2:	
					for(i = 0; i < iteraciones; i++)
					{
						total_i=omp_get_wtime();
						gate_por_canal(f_pulsos, canal_v, canal_h);

						autocorrelacion(canal_v, canal_h, autoc_v, autoc_h);

						/*Se guardan los resultados de la autocorrelacion en un archivo binario*/
						f_out=fopen(fout_name,"wb");
						fwrite(&cant_gates,1,sizeof(int),f_out); /* Indico cantidd de valores contenidosen archivo */
						fwrite(autoc_v,500,sizeof(float),f_out); /* Guardo directamente todos los resultados de canal_v */
						fwrite(autoc_h,500,sizeof(float),f_out); /* Guardo directamente todos los resultados de canal_h */
						/*fwrite(&cant_gates,1,sizeof(int),f_out); Indico cantidd de valores contenidosen archivo */
						/*for(i=0; i< GATES; i++){*/
							/*fwrite(&(autoc_v[i]),1,sizeof(float),f_out);  Guardo directamente todos los resultados de canal_v */
							/*fwrite(&(autoc_h[i]),1,sizeof(float),f_out);  Guardo directamente todos los resultados de canal_h */
						/*}*/
						fclose(f_out);
						total_f=omp_get_wtime();
						t_est = t_est + (total_f-total_i);
						im++;
						printf("Tiempo de demora total de iteracion %d: %.5f \n",i+1, total_f-total_i);
					}
					printf("\nArchivo binario de salida con resultados %s\n", fout_name);
					printf("\n Ultimo tiempo de demora total: %.5f \n", total_f-total_i);
					break;

			case 3: ver_estadisticas(t_est, im);
                	break;

			case 4: abrir_binario_resultados(fout_name);
					break;

			case 5: /*ver_resultados(autoc_v, autoc_h);*/
					for(k = 0; k < GATES; k++) {
				    	printf("R_V[%d] = %.10f, R_H[%d] = %.10f\n", k, autoc_v[k], k, autoc_h[k]);
				  	}
					break;

			case 6: break;
		}
	}

  	return 0;
}


/**/
void gate_por_canal(char *f_name, float canal_v[GATES][PULSOS], float canal_h[GATES][PULSOS]){
	FILE *f_in;
	int i,data_sample, algoritmo=1;
	int gate=0, pulso=0;
	int muestras_gate, muestras_gate_totali;
	float acumulador_v;
	float acumulador_h;
	uint16_t samples;
	float muestras_gate_totalf,error_por_muestra, error_total;
	float *cv;
	float *ch;
	float *valores;


	if((f_in= fopen(f_name,"rb")) == NULL)  /*se abre archivo en modo lectura */
	{
		perror("Fallo apertura de fichero");
		exit(1);
	}
  
  while(fread(&samples,1,sizeof(uint16_t), f_in)){
    /*La cantidad total de datos por muestra de cada pulso es 4 veces el valor samples.*/

    data_sample=4*samples;
    /*Reservo la cantidad necesaria para guardar un pulso.*/
    valores=malloc(data_sample*sizeof(float));
		cv=malloc(samples*sizeof(float));
		ch=malloc(samples*sizeof(float));
    fread(valores, data_sample, sizeof(float), f_in);
    /*indice para incrementar de a 1 la variable donde guarda los valores complejos*/
    for (i = 0; i < (2*samples); i+=2){
      /*Guarda los valores absolutos del complejo de cada canal sobre vectores diferentes.*/
			cv[i/2]=sqrt(pow(*(valores+i),2)+pow(*(valores+i+1),2));
			ch[i/2]=sqrt(pow(*(valores+i+(2*samples)),2)+pow(*(valores+i+(2*samples)+1),2));
    }

    acumulador_v=0;
    acumulador_h=0;
    muestras_gate=0;
    muestras_gate_totalf= (0.5*samples)/250; /* muestras de gate total original */
    muestras_gate_totali= (int)muestras_gate_totalf; /* muestras de gate total cortado a entero (no redondea) */
    error_por_muestra = muestras_gate_totalf - (float)muestras_gate_totali;
    
    if(error_por_muestra < 0.5) /* Redondeo hacia abajo */
    	algoritmo = 1;
    
    else if(error_por_muestra == 0) /* No hace falta redondeo */
    	algoritmo = 2;
    
    else{ /* Redondeo hacia arriba */
    	muestras_gate_totali = muestras_gate_totali + 1;
    	error_por_muestra = (float)muestras_gate_totali - muestras_gate_totalf;
    	algoritmo = 3;
    }
    
    gate=0;
    error_total = 0;

    switch(algoritmo)
    {
    	case 1:	
		    	for (i = 0; i < samples; i++) {
					acumulador_v= acumulador_v + cv[i];
					acumulador_h= acumulador_h + ch[i];
					muestras_gate++;
					if(muestras_gate==muestras_gate_totali){
						if(gate<GATES){
							error_total = error_total + error_por_muestra;
							if(error_total >= 1)
							{
								error_total = error_total - 1;
								/* Corregir error considerando una muestra mas en esta vuelta*/
								i=i+1;
								acumulador_v= acumulador_v + cv[i];
								acumulador_h= acumulador_h + ch[i];
								canal_v[gate][pulso]=acumulador_v/(muestras_gate_totali+1);
								canal_h[gate][pulso]=acumulador_h/(muestras_gate_totali+1);
							}
							else
							{
								canal_v[gate][pulso]=acumulador_v/muestras_gate_totali;
								canal_h[gate][pulso]=acumulador_h/muestras_gate_totali;
							}
							gate++;
							acumulador_v=0;
							acumulador_h=0;
							muestras_gate=0;
						}
					}
				}

    			break;
    	case 2:
				for (i = 0; i < samples; i++) {
					acumulador_v= acumulador_v + cv[i];
					acumulador_h= acumulador_h + ch[i];
					muestras_gate++;
					if(muestras_gate==muestras_gate_totali){
						if(gate<GATES){
							canal_v[gate][pulso]=acumulador_v/muestras_gate_totali;
							canal_h[gate][pulso]=acumulador_h/muestras_gate_totali;
		          gate++;
		          acumulador_v=0;
		          acumulador_h=0;
		          muestras_gate=0;
						}
					}
				}
				break;
    	case 3:
    			for (i = 0; i < samples; i++) {
					acumulador_v= acumulador_v + cv[i];
					acumulador_h= acumulador_h + ch[i];
					muestras_gate++;
					if(muestras_gate==muestras_gate_totali){
						if(gate<GATES){
							error_total = error_total + error_por_muestra;
							if(error_total >= 1)
							{
								error_total = error_total - 1;
								/* Corregir error considerando una muestra mas en esta vuelta*/
								acumulador_v= acumulador_v - cv[i];
								acumulador_h= acumulador_h - ch[i];
								canal_v[gate][pulso]=acumulador_v/(muestras_gate_totali-1);
								canal_h[gate][pulso]=acumulador_h/(muestras_gate_totali-1);
								i=i-1;
							}
							else
							{
								canal_v[gate][pulso]=acumulador_v/muestras_gate_totali;
								canal_h[gate][pulso]=acumulador_h/muestras_gate_totali;
							}
							gate++;
							acumulador_v=0;
							acumulador_h=0;
							muestras_gate=0;
						}
					}
				}
    			break;
    }
	pulso=pulso+1;

  }
}

/*Realiza una operación de autocorrelacion para las matriz v y la matriz h*/
void autocorrelacion(float canal_v[GATES][PULSOS], float canal_h[GATES][PULSOS], float p_autoc_v[], float p_autoc_h[]){
	int i=0,j=0;
	float temp_v=0;
	float temp_h=0;
	float M = 800;
	double autoc_i, autoc_f;
	autoc_i=omp_get_wtime();
	for (i = 0; i < GATES; i++) {
		for (j = 0; j < PULSOS-1; j++) {
			temp_v= temp_v + (canal_v[i][j])*(canal_v[i][j+1]);
			temp_h= temp_h + (canal_h[i][j])*(canal_h[i][j+1]);
		}
		p_autoc_v[i]=temp_v/M;
		p_autoc_h[i]=temp_h/M;
		temp_v=0;
		temp_h=0;
	}
	autoc_f=omp_get_wtime();
	/*printf("\nValores de autoc H: %f \n",p_autoc_h[0]);
	printf("\nValores de autoc V: %f \n",p_autoc_v[0]);*/
	printf("\nTiempo autocorrelación: %.5f\n", autoc_f-autoc_i);
}

/* Abre archivo con resultados y los muestra de manera adecuada*/
int abrir_binario_resultados(char *f_name)
{
	float *datos;
	int i=0, cantidad=1;
	FILE *f_in;

	f_in=fopen(f_name,"rb");

	while(fread(&cantidad,1,sizeof(int),f_in))
	{
		printf("cantidad %d\n", cantidad);
		datos=malloc(2*GATES*sizeof(float));
	    fread(datos, 2*GATES, sizeof(float), f_in);
	}
	for (i = 0; i < GATES; i++) {
    	printf("R_V[%d] = %.10f, R_H[%d] = %.10f\n", i, datos[i], i, datos[GATES+i]);
    	/*printf("cantidad %d, valores v %.10f, valores h %.10f\n", algo, *(otro+i), otro[GATES+i]);*/
  	}
	fclose(f_in);
	return 1;
}

/* Muestra resultados de la autocorrelacion sobre los diferentes canales enviados por parametro*/
void ver_resultados(float autoc_v[], float autoc_h[])
{
	int i;
	for (i = 0; i < GATES; i++) {
    	printf("R_V[%d] = %.10f, R_H[%d] = %.10f\n", i, autoc_v[i], i, autoc_h[i]);
  	}
}

/**/
/*
* Descripcion:
* @param .
*/
void ver_estadisticas(double tp_est, int im)
{
  int j =0, fa1; /* fa: flag algoritmo 1*/;
  fa1=0;

  printf("\n-----------------------------------------------------------------------------");
  printf( "\n   >>>_              ESTADISTICAS DE ALGORITMO               _<<<\n");
  printf("-----------------------------------------------------------------------------\n");

  /* METODO 1*/
	if(tp_est>0 && im>0)
	{
	  if(j==0)
	  {
	    printf("\n-----------------------------------------------------------------------------");
	    printf( "\n   >>>_                      ALGORITMO 1                      _<<<\n");
	    printf("-----------------------------------------------------------------------------\n\n");
	    printf("      |  Cores  |  Tiempo  |  Muestras  |\n");
	    j++; fa1++;
	  }
	  printf("      -------------------------------------------------------------\n");
	  printf("      |  %*.d    |  %0.4lf  |  %5.d     |\n",3, 1, (tp_est/im), im);
	}

  if(fa1==0){
    printf("\n-----------------------------------------------------------------------------");
    printf("\n   > Sin estadisticas para el algoritmo.");}

  printf("\n-----------------------------------------------------------------------------\n\n");
}

/*
* Descripcion:
* @param char *pc: puntero a l primera posicion de una cadena.
*/
int convertir_cadena(char *pc)
{
  int num_ent=1; /* Numero entero*/
  num_ent = atoi(pc);
  //printf("Numero entero: %i\n", num_ent);
  if(num_ent<1)
    num_ent=1;

  return num_ent;
}

/* time: invesigar profiling*/