#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include <mpi.h>

int numerico(char **argv){
    for(int j=1; j<3; j++){
        int i=0;
        while (argv[j][i] != '\0') {
            if (!isdigit(argv[j][i])) {
                return 1;
            }
            i++;
        }
    }
    return 0;
}

void main(int argc,char *argv[]){
    
    unsigned char sin_fichero=1;
    int nproces, myrank, k, i, m, N, j, resto, div, *cant_filas, *pos_fila;
    MPI_Status status;
    double tiempo, inicio, fin, euclides, euclides_aux, suma, normalizado, normalizado_aux;
    FILE *f;
    char *fichero;
    double **matriz, **matriz_aux, *x, *x_aux, *x_2;
    
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&nproces);
    MPI_Comm_rank(MPI_COMM_WORLD,&myrank);

    //asignar valores a número de iteraciones y tamaño de la matriz
    m=atoi(argv[1]);
    N=atoi(argv[2]);

    //inicio de la cuenta para calcular duración del programa
    //inicio = MPI_Wtime();

    if(myrank==0){
        //comprobación de parámetros
        if(argc != 4 && argc != 3){
            printf("Es necesario poner como argumento:\n 1- numero de iteraciones(numero entero)\n 2- tamano de la matriz(numero entero)\n 4- fichero\n");
            MPI_Finalize();
            return;
        }
        if(numerico(argv) == 1){
            printf("Alguno de los valores que deberia ser numerico no lo es");
            MPI_Finalize();
            return;
        }

        //asignarle espacio a la matriz original
        matriz = (double**)malloc(sizeof(double*)*N);
        for(int i = 0; i < N; i++){
            matriz[i] = (double*)malloc(sizeof(double)*N);
        }

        //comprobación de existencia de fichero
        if(argc==4){
            fichero=argv[3];
            //si el archivo no existe se crea la matriz aleatoriamente
            f=fopen(fichero,"rb");
            if (f == NULL) {
                printf("El archivo no existe, se generara una matriz aleatoriamente\n");
                sin_fichero=1;
            }else{
                fseek(f, 0, SEEK_END);
                long tam_fichero = ftell(f)/sizeof(double);
                //si el fichero seleccionado tiene menos datos de los que necesita la matriz, esta se crea aleatoriamente
                if(tam_fichero < N*N){
                    printf("No hay suficientes datos en el archivo seleccionado para generar la matriz, se generara una matriz aleatoriamente\n");
                    sin_fichero=1;
                }else{
                    //se copian los datos del fichero en la matriz
                    sin_fichero=0;
                    fseek(f, 0, SEEK_SET);
                    for(i = 0; i < N; i++){
                        fread(matriz[i],sizeof(double),N,f);
                    }
                }
                fclose(f);
            }
        }
        srand(time(NULL));
        //se genera la matriz aleatoriamente si no hay fichero viable
        if(sin_fichero==1){
            for(i = 0; i < N; i++){
                for(j = 0; j < N; j++){
                    if(i==j){
                        matriz[i][j] = 1;
                    }else{
                        matriz[i][j] =  ((double)rand() / (double)RAND_MAX) * 0.02 - 0.01;
                    }
                }
            }
        }

    }
    /*
    MPI_Barrier( MPI_COMM_WORLD);
    //inicio de la cuenta para calcular duración del programa
    inicio = MPI_Wtime();
    */
    //calcular cantidad de filas por proceso y la posición en la matriz de la primera fila de cada proceso
    div = N/nproces;
    resto = N%nproces;
    cant_filas = (int*)malloc(nproces * sizeof(int));
    pos_fila = (int*)malloc(nproces * sizeof(int));
    i=0;

    while(i<nproces){
        
        if(i<resto){
            cant_filas[i]=div+1;
        }else{
            cant_filas[i]=div;
        }
        if(i==0){
            pos_fila[i]=0;
        }else{
            pos_fila[i]=pos_fila[i-1]+cant_filas[i-1];
        }
        i++;
    }
    

    if(myrank > 0){
        matriz = (double**)malloc(sizeof(double*)*cant_filas[myrank]);
        for(i=0; i<cant_filas[myrank]; i++){
            matriz[i] = (double*)malloc(sizeof(double)*N);
        }
    }
    //asignarle espacio al vector x e igualarlo a 1.0
    x = (double*)malloc(sizeof(double)*N);
    x_aux = (double*)malloc(sizeof(double)*cant_filas[myrank]);
    x_2 = (double*)malloc(sizeof(double)*cant_filas[myrank]);
    for(i=0; i<N; i++){
        x[i]=1.0;
    }

    if(nproces>1){
        if(myrank==0){
            for(i = 1; i<nproces; i++){
                k=pos_fila[i];
                for(j = 0; j<cant_filas[i]; j++){
                    MPI_Send(matriz[k],N,MPI_DOUBLE,i,1,MPI_COMM_WORLD);
                    k++;
                }
            }
        }else{
            for(j = 0; j<cant_filas[myrank]; j++){
                MPI_Recv(matriz[j],N,MPI_DOUBLE,0,1,MPI_COMM_WORLD, &status);
                
            }
            
        }
    }
    
    //iteraciones para realizar los calculos necesarios
    for(i=0; i<m; i++){

        //se multiplica el vector por la matriz y se le resta el resultado al vector resultante de la iteración anterior
        for(k = 0; k < cant_filas[myrank]; k++){
            suma=0;

            for(j = 0; j < N; j++){
                suma+=matriz[k][j]*x[j];
            }

            x_aux[k] = suma;
            if(i>0){
                x_aux[k] = x[pos_fila[myrank]+k]-x_aux[k];
            }
        
        }
        //realizamos euclides
        euclides=0;
        for(j=0; j<cant_filas[myrank]; j++){
            euclides+=pow(x_aux[j],2);
        }

        MPI_Reduce(&euclides, &euclides_aux, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if(myrank==0){
            normalizado=sqrt(euclides_aux);
            printf("Norma de iteracion %d:\n",i+1);
            printf("%lf\n", normalizado);  
        }

         MPI_Bcast(&normalizado, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        if(i<m-1){
            //terminamos de normalizar el vector
            for(j=0; j<cant_filas[myrank]; j++){
                x_2[j]=x_aux[j]/normalizado; 
            }  

            MPI_Allgatherv(x_2, cant_filas[myrank], MPI_DOUBLE, x, cant_filas, pos_fila, MPI_DOUBLE, MPI_COMM_WORLD);
        }
    }    
    
    //liberar memoria de los vectores
    free(x_aux);
    free(x);
    free(x_2);

    if(myrank==0){
        
        //liberar espacio de la matriz
        for(i = 0; i<N; i++){
            free(matriz[i]);
        }
        free(matriz);
        
    }
    /*
    //fin de medición de tiempo
    fin = MPI_Wtime();
    tiempo = (fin - inicio);
    
    printf("\nTiempo medido: %.3lf segundos del proceso %d.\n", tiempo, myrank);
    */

    MPI_Finalize();
    return;
}