#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <string.h>
#include <unistd.h>

#define GATES 500
#define PULSOS 800 /** 100 pulsos por 8 grado de acimut considerados (100*8=800)*/
#define name_size 100
#define size_comando 256
#define MAX_HILOS 100 /** Soporte maximo de hilos para analisis de estadisticas*/

/** Declaracion de funciones utilizadas en aplicacion */
int abrir_binario_resultados(char *f_name);
int convertir_cadena(char *pc);
void ver_estadisticas_procedural(float tp_est, int im);
void ver_estadisticas_paralelo(float tp_est[], int im[], float speedUp[], float porcentajes[]);
void borrar_estadisticas(float *t_est, float tp_est[], int *imp, int im_parallel[], float speedup[], float porcentajes[]);


/** Funcion principal */
int main(int argc, char *argv[]) {

	/** Nombres de ficheros utilizados*/
	char fout_name[name_size]="resultados_adpc"; /* resultados de autocorrelacion discreta por pulsos de cada canal*/
	char f_pulsos[name_size]="pulsos.iq";

	/** Variables para contar tiempo de ejecución*/
	float total_i, total_f;

	/** Varibles asociadas a parametros generales en algoritmo procedural y paralelo*/
	int opcion =0;
	int i;
	int iteraciones=1;
	char cadena[size_comando];
	char analisis_cadena[size_comando]="-";

	/** Variables vinculadas a parametros de algoritmo procedural*/
	char comando_procedural_modo1[size_comando]="./TP2_procedural.out";
	char comando_procedural_modo2[size_comando]="./TP2_procedural.out 2"; /** con este comando se muestran los resultados de autocorrelacion*/
	float t_est=0; /** tiempo estadistico segun numero de iteraciones por algoritmo procedural*/
	int imp=0; /** iterador de metodo procedural. Es im_parallel[0] en el vector im_parallel*/

	/** Variables vinculadas a parametros de algoritmo paralelao*/
	char comando_paralelo[size_comando]="./TP2_paralelo_2.out";
	char comando_paralelo_modo1[size_comando]="./TP2_paralelo_2.out 2";
	char comando_paralelo_modo2[size_comando]="./TP2_paralelo_2.out 2 1";
	int hilos=1;
	float tp_estaux=0;
	int im_paux=0;
	int im_parallel[MAX_HILOS]={0}; /** Iterador de metodo paraleloslo, 1 por cada numero de hilo utilizado, soporte a MAX_HILOS*/
	float tp_est[MAX_HILOS]={0};
	float speedUp[MAX_HILOS]={0}, porcentajes[MAX_HILOS]={0};

	while(opcion!=11)
	{
		printf( "\n   >>>_			MENU 			_<<<\n");
        printf( "\n   1. Cambiar numero de iteraciones por algoritmo.");
        printf( "\n   2. Cambiar numero de hilos paralelos.");
        printf( "\n   3. Ejecutar algoritmo procedural.");
        printf( "\n   4. Ejecutar algoritmo procedural con solicitud de resultados.");
        printf( "\n   5. Ejecutar algoritmo paralelo.");
        printf( "\n   6. Ejecutar algoritmo paralelo con solicitud de resultados.");
        printf( "\n   7. Ver estadisticas de algoritmo procedural.");
        printf( "\n   8. Ver estadisticas de algoritmo paralelo.");
        printf( "\n   9. Borrar estadisticas.");
        printf( "\n   10. Mostrar resultados desde archivo binario %s.", fout_name);
        printf( "\n   11. Salir." );
        printf( "\n\n   Introduzca opcion (1-11): ");

        scanf( "%d", &opcion );
		switch(opcion)
		{
			case 1:	/** Cambiar numero de iteraciones por algoritmo. */
					while(strcmp(analisis_cadena,"-")==0)
	        		{
						memset(cadena, '\0', 256);
						memset(analisis_cadena, '\0', size_comando);
		                printf( "\n   > Introduzca el numero de iteraciones deseado: ");
		                scanf( "%s", cadena);
		                strncat(analisis_cadena, cadena, 1);
		                if(strcmp(analisis_cadena,"-")==0)
		                	printf("\t> Nose permiten numeros negativos intente nuevamente\n");
		            }
	                iteraciones = convertir_cadena(cadena);
	                printf("\t> Se establecieron %d iteraciones.\n",iteraciones );
	                memset(analisis_cadena, '\0', size_comando);
					strcpy(analisis_cadena, "-"); /* Reestablesco valor de analisis*/
	                break;

	        case 2:	/** Cambiar numero de hilos paralelos. */
	        		while(strcmp(analisis_cadena,"-")==0)
	        		{
		        		memset(cadena, '\0', size_comando);
		        		memset(analisis_cadena, '\0', size_comando);
		                printf( "\n   > Introduzca el numero de hilos paralelos deseados: ");
		                scanf( "%s", cadena);
		                strncat(analisis_cadena, cadena, 1);
		                if(strcmp(analisis_cadena,"-")==0)
		                	printf("\t> No se permiten numeros negativos intente nuevamente\n");
		            }
	                hilos = convertir_cadena(cadena);
					memset( comando_paralelo_modo1, '\0', size_comando);
					memset( comando_paralelo_modo2, '\0', size_comando);
					strcat(comando_paralelo_modo1, comando_paralelo);
					strcat(comando_paralelo_modo1, " ");
					strcat(comando_paralelo_modo1, cadena);
					strcat(comando_paralelo_modo2, comando_paralelo_modo1);
					strcat(comando_paralelo_modo2, " ");
					strcat(comando_paralelo_modo2, "1");
					printf("\t* Comando 1:%s\n\t* Comando 2:%s\n\t* hilos: %d\n", comando_paralelo_modo1, comando_paralelo_modo2, hilos);
					memset(analisis_cadena, '\0', size_comando);
					strcpy(analisis_cadena, "-"); /** Reestablesco valor de analisis*/
					/** Ejemplos de comando a ejecutar en modo paralelo
					* ./TP2_paralelo-out 2 : este comando significa que se ejecuta algoritmo paralelo con 2 hilos.
					* ./TP2_paralelo-out 2 1: este es igual al anterior pero con el otro parametro indica mostrar los resultados de la autocorrelacion obtenida.
					*/
	                break;

			case 3:	/** Ejecutar algoritmo procedural. */
					printf("\n\t* El numero de pulsos que tiene el archivo %s es: %d \n",f_pulsos, PULSOS);
					printf("\t* Las componetes de rango asociadas a cada gate son: %d \n", GATES);
					printf("\t* Se utlizan estructuras de datos del tipo matriz[GATES][PULSO] = canal[%d][%d] \n\n", GATES, PULSOS);
					for(i = 0; i < iteraciones; i++)
					{
						total_i=omp_get_wtime();
						/** Ejecucion de proceso con algoritmo procedural*/
						system(comando_procedural_modo1);
						total_f=omp_get_wtime();
						t_est = t_est + (total_f-total_i);
						imp++;
						printf("\t* Tiempo de demora total de iteracion %d: %.5f \n\n",i+1, total_f-total_i);
					}
					im_parallel[0] = imp;
					tp_est[0]= t_est;
					printf("\t* Archivo binario de salida con resultados %s\n\n", fout_name);
					break;

			case 4:	/** Ejecutar algoritmo procedural con solicitud de resultados. */
					total_i=omp_get_wtime();
					/** Ejecucion de proceso con algoritmo procedural*/
					system(comando_procedural_modo2);
					total_f=omp_get_wtime();
					t_est = t_est + (total_f-total_i);
					imp++;
					im_parallel[0] = imp;
					tp_est[0]= t_est;
					printf("Tiempo de demora total de iteracion %d: %.5f \n",i+1, total_f-total_i);
					printf("\nArchivo binario de salida con resultados %s\n", fout_name);
					break;

			case 5:	/** Ejecutar algoritmo paralelo */
					printf("\n\t* El numero de pulsos que tiene el archivo %s es: %d \n",f_pulsos, PULSOS);
					printf("\t* Las componetes de rango asociadas a cada gate son: %d \n", GATES);
					printf("\t* Se utlizan estructuras de datos del tipo matriz[GATES][PULSO] = canal[%d][%d] \n\n", GATES, PULSOS);
					tp_estaux=0;
					im_paux=0;
					for(i = 0; i < iteraciones; i++)
					{
						total_i=omp_get_wtime();
						/** Ejecucion de proceso con algoritmo procedural*/
						system(comando_paralelo_modo1);
						total_f=omp_get_wtime();
						tp_estaux = tp_estaux + (total_f-total_i);
						im_paux++;
						printf("Tiempo de demora total de iteracion %d: %.5f \n\n",i+1, total_f-total_i);
					}
					tp_est[hilos] = tp_est[hilos] + tp_estaux;
					im_parallel[hilos]= im_parallel[hilos] + im_paux;
					printf("\nArchivo binario de salida con resultados %s\n\n", fout_name);
					break;

			case 6:	/** Ejecutar algoritmo paralelo con solicitud de resultados. */
					total_i=omp_get_wtime();
					/** Ejecucion de proceso con algoritmo procedural*/
					system(comando_paralelo_modo2);
					total_f=omp_get_wtime();
					tp_est[hilos] = tp_est[hilos] + (total_f-total_i);
					im_parallel[hilos]++;
					printf("Tiempo de demora total de iteracion %d: %.5f \n",i+1, total_f-total_i);
					printf("\nArchivo binario de salida con resultados %s\n", fout_name);
					break;


			case 7: /** Ver estadisticas de algoritmo procedural. */
					ver_estadisticas_procedural(t_est, imp);
                	break;

            case 8: /** Ver estadisticas de algoritmo paralelo. */
            		ver_estadisticas_paralelo(tp_est, im_parallel, speedUp, porcentajes);
                	break;

            case 9: /** Borrar estadisticas. */
                	printf("\t> Se eliminaran todos los datos asociados a las estadistica. Continuar(si/no)?\n");
                	memset(cadena, '\0', 256);
		            scanf( "%s", cadena);
		            if(strcmp(cadena, "si")==0)
		            {
	                	borrar_estadisticas(&t_est, tp_est, &imp, im_parallel, speedUp, porcentajes);
	                	printf("\t> Se han eliminado los datos de estadisticas\n\n");
	                }
					break;

			case 10: /** Mostrar resultados desde archivo binario. */ 
					abrir_binario_resultados(fout_name);
					break;

			case 11: /** Salir. */ 
					break;
		}
	}

  	return 0;
}


/** Abre archivo con resultados y los muestra de manera adecuada*/
int abrir_binario_resultados(char *f_name)
{
	double *datos;
	int i=0, cantidad=1;
	FILE *f_in;

	f_in=fopen(f_name,"rb");

	while(fread(&cantidad,1,sizeof(int),f_in))
	{
		printf("cantidad %d\n", cantidad);
		datos=malloc(2*GATES*sizeof(double));
	    	fread(datos, 2*GATES, sizeof(double), f_in);
	}
	for (i = 0; i < GATES; i++) {
    		printf("R_V[%d] = %.10f, R_H[%d] = %.10f\n", i, datos[i], i, datos[GATES+i]);
  	}
	fclose(f_in);
	free(datos);
	return 1;
}


/**
* Descripcion: Muestra las estadisticas del proceso en procedural
*/
void ver_estadisticas_procedural(float tp_est, int im)
{
  int j =0, fa1; /** fa: flag algoritmo 1*/;
  fa1=0;

  printf("\n-----------------------------------------------------------------------------");
  printf( "\n   >>>_              ESTADISTICAS DE ALGORITMO               _<<<\n");
  printf("-----------------------------------------------------------------------------\n");

  /** METODO 1 */
  if(tp_est>0 && im>0)
   {
	  if(j==0)
	  {
	    printf("\n-----------------------------------------------------------------------------");
	    printf( "\n   >>>_                      ALGORITMO PROCEDURAL                      _<<<\n");
	    printf("-----------------------------------------------------------------------------\n\n");
	    printf("      |  Cores  |  Tiempo  |  Muestras  |\n");
	    j++; fa1++;
	  }
	  printf("      -------------------------------------------------------------\n");
	  printf("      |  %*.d    |  %0.4f  |  %5.d     |\n",3, 1, (tp_est/im), im);
    }

  if(fa1==0){
    printf("\n-----------------------------------------------------------------------------");
    printf("\n   > Sin estadisticas para el algoritmo procedural.");}

  printf("\n-----------------------------------------------------------------------------\n\n");
}


/**
* Descripcion:  Muestra las estadisticas en paralelo.
*/
void ver_estadisticas_paralelo(float tp_est[], int im[], float speedUp[], float porcentajes[])
{
	int i, j =0, fa1; /** fa: flag algoritmo 1*/
	char sim=37; /** simbolo de caracter equvalente a %*/
	fa1=0;

	printf("\n-----------------------------------------------------------------------------");
	printf( "\n   >>>_              ESTADISTICAS DE ALGORITMOS               _<<<\n");
	printf("-----------------------------------------------------------------------------\n");

	/** METODO 1. */
	for(i=0; i<MAX_HILOS-1; i++)
	{
		if(tp_est[i]>0 && im[i]>0)
		{
	  		if(j==0)
	  		{
	    			printf("\n-----------------------------------------------------------------------------");
	    			printf( "\n   >>>_                      ALGORITMO PARALELO                      _<<<\n");
			        printf("-----------------------------------------------------------------------------\n\n");
		                printf("      |  Cores  |  Tiempo  | SpeedUp  |  Porcentaje  |  Muestras  |\n");
		          	j++;
				fa1++;
	  		}
	  	/** Calculo de speedUp practico*/
	  	speedUp[i] = (tp_est[0]/im[0])/(tp_est[i]/im[i]);
	  	porcentajes[i] = (speedUp[i]*100)-100;
	  	printf("      -------------------------------------------------------------\n");
	  	if(i==0)
	  		printf("      |  %*.d    |  %0.4f  |  %0.4f  |  %5.0f%c      |  %5.d     |\n",3, i+1, (tp_est[i]/im[i]), speedUp[i], porcentajes[i], sim, im[i]);
	  
	  	else
	  		printf("      |  %*.d    |  %0.4f  |  %0.4f  |  %5.0f%c      |  %5.d     |\n",3, i, (tp_est[i]/im[i]), speedUp[i], porcentajes[i], sim, im[i]);
		}
	}

	if(fa1==0){
	printf("\n-----------------------------------------------------------------------------");
	printf("\n   > Sin estadisticas para el algoritmo paralelo.");}

	printf("\n-----------------------------------------------------------------------------\n\n");
}

/**
* Descripcion:
* @param char *pc: puntero a l primera posicion de una cadena.
*/
int convertir_cadena(char *pc)
{
	int num_ent=1; /** Numero entero*/
	num_ent = atoi(pc);
	if(num_ent<1){
		num_ent=1;
		memset(pc, '\0', 256);
		strcpy(pc, "1");
	}
	
	return num_ent;
}

void borrar_estadisticas(float *t_est, float tp_est[], int *imp, int im_parallel[], float speedup[], float porcentajes[])
{
	int i=0;
	*t_est=0;
	*imp=0;
	for(i=0; i<MAX_HILOS-1; i++)
	{
		tp_est[i]=0;
		im_parallel[i]=0;
		speedup[i]=0;
		porcentajes[0]=0;
	}

}


