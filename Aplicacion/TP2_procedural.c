#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <unistd.h>
#include <omp.h>
#include <string.h>

#define GATES 500
#define PULSOS 100*8 /* 100 pulsos por 8 grado de acimut considerados (100*8=800)*/
#define name_size 100

/* Declaracion de funciones utilizadas en aplicacion */
void gate_por_canal(char *f_name, float canal_v[GATES][PULSOS], float canal_h[GATES][PULSOS]);
void autocorrelacion(float canal_v[GATES][PULSOS], float canal_h[GATES][PULSOS], double p_autoc_v[], double p_autoc_h[]);

/* Funcion principal */
int main(int argc, char const *argv[]) {

	/* Nombres de ficheros utilizados*/
	char f_pulsos[name_size]="pulsos.iq";
	char fout_name[name_size]="resultados_adpc"; /* resultados de autocorrelacion discreta por pulsos de cada canal*/

	
	/*Variables para contar tiempo de ejecución*/
	double total_i, total_f;

	/* Matrices de variables para las estructura gate-pulso por canal*/
	float canal_v[GATES][PULSOS]; /* inicializada a cero en todos sus valores*/
	float canal_h[GATES][PULSOS]; /* inicializada a cero en todos sus valores*/
	int i, j;
	/* Inicializacion de matriz*/
	for(i=0; i<GATES; i++)
		for(j=0; j<PULSOS; j++)
			canal_v[i][j]=canal_h[i][j]=0;
	
	/* Vectores con los resultados de la autocorrelacion discreta por pulsos de cada canal*/
	double autoc_v[GATES]={0};
	double autoc_h[GATES]={0};

	int k=0;
	
	/* Variables de archivos.*/
	FILE *f_out;
	int cant_gates= (int)GATES;

	/* Inicio toma de tiempos por calculos de autocorrelacion*/
	total_i=omp_get_wtime();
	gate_por_canal(f_pulsos, canal_v, canal_h);

	autocorrelacion(canal_v, canal_h, autoc_v, autoc_h);

	/*Se guardan los resultados de la autocorrelacion en un archivo binario*/
	f_out=fopen(fout_name,"wb");
	fwrite(&cant_gates,1,sizeof(int),f_out); /* Indico cantidd de valores contenidos en archivo */
	fwrite(autoc_v,500,sizeof(double),f_out); /* Guardo directamente todos los resultados de canal_v */
	fwrite(autoc_h,500,sizeof(double),f_out); /* Guardo directamente todos los resultados de canal_h */
	/*fwrite(&cant_gates,1,sizeof(int),f_out); Indico cantidd de valores contenidosen archivo */
	/*for(i=0; i< GATES; i++){*/
		/*fwrite(&(autoc_v[i]),1,sizeof(float),f_out);  Guardo directamente todos los resultados de canal_v */
		/*fwrite(&(autoc_h[i]),1,sizeof(float),f_out);  Guardo directamente todos los resultados de canal_h */
	/*}*/
	fclose(f_out);
	total_f=omp_get_wtime();
	printf("\t* Tiempo requerido para calculos (calculo de GATES y autocorrelacion): %.5f \n", total_f-total_i);

	if(argc==2)
	{
		printf("\n \t* Los resultados capturados en las estructuras de datos son: \n");
		for(k = 0; k < GATES; k++) {
	    	printf("R_V[%d] = %.10g, R_H[%d] = %.10g\n", k, autoc_v[k], k, autoc_h[k]);
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

    data_sample=4*(int)samples;
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
								if(i<muestras_gate_totali)/* corregir*/
								{
									acumulador_v= acumulador_v + cv[i];
									acumulador_h= acumulador_h + ch[i];
									canal_v[gate][pulso]=acumulador_v/(muestras_gate_totali+1);
									canal_h[gate][pulso]=acumulador_h/(muestras_gate_totali+1);
								}
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
    free(cv); /* Libero espacios utilizados temporalmente, para liberar RAM*/
    free(ch);
    free(valores);
	pulso=pulso+1;

  }
  fclose(f_in);
  
  /*for(i=0; i<PULSOS; i++)
  {
  	printf("Gate[%d]=%f\n", i,canal_h[499][i]);
  }*/
}

/*Realiza una operación de autocorrelacion para las matriz v y la matriz h*/
void autocorrelacion(float canal_v[GATES][PULSOS], float canal_h[GATES][PULSOS], double p_autoc_v[], double p_autoc_h[]){
	int i=0,j=0;
	double temp_v=0;
	double temp_h=0;
	double M = 800;
	double autoc_i, autoc_f;
	autoc_i=omp_get_wtime();
	for (i = 0; i < GATES; i++) {
		for (j = 0; j < PULSOS-1; j++) {
			temp_v= temp_v + (double)((canal_v[i][j])*(canal_v[i][j+1]));
			temp_h= temp_h + (double)((canal_h[i][j])*(canal_h[i][j+1]));
		}
		p_autoc_v[i]=temp_v/M;
		p_autoc_h[i]=temp_h/M;
		temp_v=0;
		temp_h=0;
	}
	autoc_f=omp_get_wtime();
	/*printf("\nValores de autoc H: %f \n",p_autoc_h[499]);
	printf("\nValores de autoc V: %f \n",p_autoc_v[0]);*/
	printf("\t* Tiempo autocorrelación: %.5f\n", autoc_f-autoc_i);
}
/* time: invesigar profiling*/