#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include <mpi.h>

//funcion para calcular el histograma
int *histograma(int filas, int columnas, unsigned char **matriz){

    int *contador;
    contador = (int*)malloc(256*sizeof(int));
    memset(contador,0,256*sizeof(int));
    //contar la cantidad de veces que se repite cada valor en la matriz
    for(int i = 1; i<filas+1; i++){
        for(int j = 1; j<columnas+1; j++){
            contador[matriz[i][j]]++;
        }
    }
    return contador;
}

void copiar_hist(int *contador){
    FILE *f;
    //abrir el archivo de texto y copiar los datos del histograma
    f = fopen("histograma.txt","w");
    if (f == NULL) {
        printf("No se puede abrir el archivo\n");
        return;
    }
    for(int i=0; i<256; i++){
        fprintf(f, "%d - %d\n", i, contador[i]);
    }
    fclose(f);
    return;
}

//función para calcular la matriz de Medias
void calculoMedia(unsigned char **matriz, int filas, int columnas, int myrank){
    int suma;
    //asignar espacio para la nueva matriz de calculo de media
    unsigned char **matrizMedia;
    matrizMedia = (unsigned char**)malloc((filas) * sizeof(unsigned char*));
    for(int i = 0; i<filas; i++){
        matrizMedia[i] = (unsigned char*)malloc((columnas) * sizeof(unsigned char));
    }
    //calcular la media de cada dato de la matriz original y guardarlos en la matriz de medias
    for(int i = 0; i<filas; i++){
        for(int j = 0; j<columnas; j++){
            suma = 0;
            for(int k = 0; k<3; k++){
                for(int l = 0; l<3; l++){
                    //realizar la suma de los 9 elementos de la minimatriz
                    suma += matriz[i+k][j+l];
                }
            }
            //dividir la suma de la minimatriz entre 9 para hayar la media
            matrizMedia[i][j] = suma / 9;        
        }
    }
    
    //enviar datos de la matriz de medias al root
    for(int i=0; i<filas; i++){
        MPI_Send(matrizMedia[i],columnas,MPI_UNSIGNED_CHAR,0,3,MPI_COMM_WORLD);
    }

    return;
}

//función para calcular la matriz de Sobel
void sobel(unsigned char **matriz, int filas, int columnas, int myrank, int nproces){
    //asignar espacio para la matriz de resultados de C y F, en este caso son ints porque hay resultados que pueden ser negativos
    int c;
    int f;
    //declarar las matrices C y F
    int C[3][3] = {{-1, 0, 1},{-2, 0, 2},{-1, 0, 1}};
    int F[3][3] = {{-1, -2, -1},{0, 0, 0},{1, 2, 1}};
    //asignar espacio para la nueva matriz de calculo de media
    unsigned char **matrizSobel;
    if(myrank ==0 && myrank == nproces-1){//si hay un solo proceso
        matrizSobel = (unsigned char**)malloc((filas-2) * sizeof(unsigned char*));
        for(int i = 0; i<filas-1; i++){
            matrizSobel[i] = (unsigned char*)malloc((columnas-2) * sizeof(unsigned char));
        }
        //realizar los calculos necesarios
        for(int i = 1; i<filas-1; i++){
            for(int j = 1; j<columnas-1; j++){
                c = 0;
                f = 0;
                for(int k = 0; k<3; k++){
                    for(int l = 0; l<3; l++){
                        c += matriz[i+k][j+l]*C[k][l];
                        f += matriz[i+k][j+l]*F[k][l];  
                    }
                }
                matrizSobel[i-1][j-1] = sqrt( c*c + f*f );
            }
        }
        
    }else if(myrank ==0 || myrank == nproces-1){//si es el proceso 0 o el último
        matrizSobel = (unsigned char**)malloc((filas-1) * sizeof(unsigned char*));
        for(int i = 0; i<filas-1; i++){
            matrizSobel[i] = (unsigned char*)malloc((columnas-2) * sizeof(unsigned char));
        }
        if(myrank ==0){
            for(int i = 1; i<filas; i++){
                for(int j = 1; j<columnas-1; j++){
                    c = 0;
                    f = 0;
                    for(int k = 0; k<3; k++){
                        for(int l = 0; l<3; l++){
                            c += matriz[i+k][j+l]*C[k][l];
                            f += matriz[i+k][j+l]*F[k][l];
                            
                        }
                    }
                    matrizSobel[i-1][j-1] = sqrt( c*c + f*f );
                }
            }
        }else{
            for(int i = 0; i<filas-1; i++){
                for(int j = 1; j<columnas-1; j++){
                    c = 0;
                    f = 0;
                    for(int k = 0; k<3; k++){
                        for(int l = 0; l<3; l++){
                            c += matriz[i+k][j+l]*C[k][l];
                            f += matriz[i+k][j+l]*F[k][l];
                            
                        }
                    }
                    matrizSobel[i][j-1] = sqrt( c*c + f*f );
                }
            }
        }
        
    }else {//resto de procesos
        matrizSobel = (unsigned char**)malloc((filas) * sizeof(unsigned char*));
        for(int i = 0; i<filas; i++){
            matrizSobel[i] = (unsigned char*)malloc((columnas-2) * sizeof(unsigned char));
        }
        //calcular la detección de bordes (sobel)  de cada dato de la matriz original y guardarlos en la matriz de medias
        for(int i = 0; i<filas; i++){
            for(int j = 1; j<columnas-1; j++){
                c = 0;
                f = 0;
                for(int k = 0; k<3; k++){
                    for(int l = 0; l<3; l++){
                        c += matriz[i+k][j+l]*C[k][l];
                        f += matriz[i+k][j+l]*F[k][l];
                        
                    }
                }
                matrizSobel[i][j-1] = sqrt( c*c + f*f );
                //printf("%d  ", matrizSobel[i-1][j-1]);   
            }
        }
    }
    
    //enviar datos de la matriz Sobel al root
    if(myrank == 0 && myrank == nproces-1){
        for(int i=0; i<filas-2; i++){
            MPI_Send(matrizSobel[i],columnas-2,MPI_UNSIGNED_CHAR,0,3,MPI_COMM_WORLD);
        }
    }else if(myrank == 0 || myrank == nproces-1){
        for(int i=0; i<filas-1; i++){
            MPI_Send(matrizSobel[i],columnas,MPI_UNSIGNED_CHAR,0,3,MPI_COMM_WORLD);
        }
    }else{
        for(int i=0; i<filas; i++){
            MPI_Send(matrizSobel[i],columnas,MPI_UNSIGNED_CHAR,0,3,MPI_COMM_WORLD);
        }
    }

    return;
}

void copiar_datos_bin(unsigned char **matriz, int filas, int columnas, char *fichero, char *tipo){
    FILE *f;
    //copiar los datos en el fichero binario
    char newFile[30]="";
    strcat(newFile,tipo);
    strcat(newFile,fichero);
    f = fopen(newFile, "wb");
    for(int i = 0; i<filas; i++){
        fwrite(matriz[i], sizeof(unsigned char), columnas, f);
    }
    fclose(f);
}

void main(int argc,char *argv[])
{
    //declaración de variables
    int nproces, myrank, k, i, err, filas, columnas, *cant_filas, *pos_fila, resto, div;
    MPI_Status status;
    double tiempo, inicio, fin;
    FILE *f;
    char *fichero;
    unsigned char procesado;
    unsigned char **matriz, **matriz_aux;
    int *contador;
    int contador_aux[256] = {0};

    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&nproces);
    MPI_Comm_rank(MPI_COMM_WORLD,&myrank);

    filas = atoi(argv[2]);
    columnas = atoi(argv[3]);
    procesado = atoi(argv[4]);

    //calcular cantidad de filas por proceso y la posición en la matriz de la primera fila de cada proceso
    div = filas/nproces;
    resto = filas%nproces;
    cant_filas = (int*)malloc(nproces * sizeof(int));
    pos_fila = (int*)malloc(nproces * sizeof(int));
    i=0;
    k=0;

    //calcular las filas que le corresponden a cada proceso y la posición de la primera fila que le corresponde
    while(i<nproces){
        if(i==nproces-1){
            cant_filas[i]=div;
        }else{
            if(i<resto){
                cant_filas[i]=div+1;
            }else{
                cant_filas[i]=div;

            }
        }
        if(i==0){
            pos_fila[i]=0;
        }else if(i==1){
            pos_fila[i]=pos_fila[i-1]+cant_filas[i-1]-1;
        }
        else{
            pos_fila[i]=pos_fila[i-1]+cant_filas[i-1];
        }
        i++;
    
    }

    //proceso root 0
    if(myrank==0){
        //comenzar a medir el tiempo que se demora
        inicio = MPI_Wtime();
        fichero = argv[1];

        //Abrimos el archivo binario .raw
        f = fopen(fichero, "rb");
        if (f == NULL) {
            printf("No se puede abrir el archivo binario\n");
            MPI_Finalize();
            return;
        }

        //Le asignamos espacio de manera dinámica a la matriz en la que guardaremos los datos
        matriz = (unsigned char**)malloc((filas) * sizeof(unsigned char*));
        for(int i = 0; i<filas; i++){
            matriz[i] = (unsigned char*)malloc((columnas) * sizeof(unsigned char));
            //copiamos los datos del fichero en la matriz
            fread(matriz[i], sizeof(unsigned char), columnas, f);
        }

        //cerramos el archivo binario
        fclose(f);

        //enviar la información de las filas de la matriz que le corresponden a cada proceso
        for(i = 0; i<nproces; i++){
            k=pos_fila[i];
            if(i==nproces-1 && i==0){
                for(int j = 0; j<cant_filas[i]; j++){
                    MPI_Send(matriz[k],columnas,MPI_UNSIGNED_CHAR,0,10,MPI_COMM_WORLD);
                    k++;
                }
            }else if(i==nproces-1 || i==0){
                for(int j = 0; j<cant_filas[i]+1; j++){
                    MPI_Send(matriz[k],columnas,MPI_UNSIGNED_CHAR,i,1,MPI_COMM_WORLD);
                    k++;
                }
            }else{
                for(int j = 0; j<cant_filas[i]+2; j++){
                    MPI_Send(matriz[k],columnas,MPI_UNSIGNED_CHAR,i,1,MPI_COMM_WORLD);
                    k++;
                }
            }
        }
        
    }
        //reservar memoria para la matriz
        matriz_aux = (unsigned char**)malloc((cant_filas[myrank]+2) * sizeof(unsigned char*));
        for(int i = 0; i<cant_filas[myrank]+2; i++){
            matriz_aux[i] = (unsigned char*)malloc((columnas+2) * sizeof(unsigned char));
        }

        //recivir la información enviada por el proceso root 0
        if(myrank == nproces-1 && myrank == 0){
            for(int j = 1; j<cant_filas[myrank]+1; j++){
                MPI_Recv(&matriz_aux[j][1],columnas,MPI_UNSIGNED_CHAR,0,10,MPI_COMM_WORLD, &status);
            }
        }else if(myrank == nproces-1){
            for(int j = 0; j<cant_filas[myrank]+1; j++){
                MPI_Recv(&matriz_aux[j][1],columnas,MPI_UNSIGNED_CHAR,0,1,MPI_COMM_WORLD, &status);
            }
        }else if(myrank == 0){
            for(int j = 1; j<cant_filas[myrank]+2; j++){
                MPI_Recv(&matriz_aux[j][1],columnas,MPI_UNSIGNED_CHAR,0,1,MPI_COMM_WORLD, &status);
            }
        }else{
            for(int j = 0; j<cant_filas[myrank]+2; j++){
                MPI_Recv(&matriz_aux[j][1],columnas,MPI_UNSIGNED_CHAR,0,1,MPI_COMM_WORLD, &status);
            }
        }
    
    //Extención simétrica: rellenar la primera y última columna para poder procesar los datos
    for(int i = 0; i<cant_filas[myrank]+2; i++){
        matriz_aux[i][0] = matriz_aux[i][2];
        matriz_aux[i][columnas+1] = matriz_aux[i][columnas-1];
    }
    if(myrank==0){
        //rellenar la primera fila para poder procesar los datos
        for(int i = 1; i<columnas+1; i++){
            matriz_aux[0][i] = matriz_aux[2][i];
        }
        matriz_aux[0][0] = matriz_aux[2][2];
        matriz_aux[0][columnas+1] = matriz_aux[2][columnas-1];
    }
    if(myrank==nproces-1){
        //rellenar la última fila para poder procesar los datos
        for(int i = 1; i<columnas+1; i++){
            matriz_aux[cant_filas[myrank]+1][i] = matriz_aux[cant_filas[myrank]-1][i];
        } 
        matriz_aux[cant_filas[myrank]+1][0] = matriz_aux[cant_filas[myrank]-1][2];
        matriz_aux[cant_filas[myrank]+1][columnas+1] = matriz_aux[cant_filas[myrank]-1][columnas-1];
    }    
    //si se pide calcular la media
    if(procesado == 1){
        calculoMedia(matriz_aux, cant_filas[myrank], columnas, myrank); 
        if(myrank==0){
            for(i = 0; i<nproces; i++){
                k=pos_fila[i];
                if(i>0){
                    k+=1;
                }
                for(int j = 0; j<cant_filas[i]; j++){
                    MPI_Recv(matriz[k],columnas,MPI_UNSIGNED_CHAR,i,3,MPI_COMM_WORLD, &status);
                    k++;
                }
            }
    
            copiar_datos_bin(matriz, filas, columnas, fichero, "Media");
        }
    }else if(procesado == 3){//si se pide calcular Sobel

        sobel(matriz_aux, cant_filas[myrank], columnas, myrank, nproces);
        if(myrank==0){
            for(i = 0; i<nproces; i++){
                k=pos_fila[i];

                //recivir la información enviada al proceso root 0
                if(i==0 && i ==nproces-1){
                    for(int j = 0; j<cant_filas[i]-2; j++){
                        MPI_Recv(&matriz[k+1][1],columnas-1,MPI_UNSIGNED_CHAR,i,3,MPI_COMM_WORLD, &status);
                        k++;
                    }
                }else if(i==0 || i ==nproces-1){
                    for(int j = 0; j<cant_filas[i]-1; j++){
                        MPI_Recv(&matriz[k+1][1],columnas,MPI_UNSIGNED_CHAR,i,3,MPI_COMM_WORLD, &status);
                        k++;
                    }
                }else{
                    for(int j = 0; j<cant_filas[i]; j++){
                        MPI_Recv(&matriz[k+1][1],columnas,MPI_UNSIGNED_CHAR,i,3,MPI_COMM_WORLD, &status);
                        k++;
                    }
                }       
            }
            //Extención simétrica: rellenar la primera y última columna para poder procesar los datos
            for(int i = 0; i<filas; i++){
                matriz[i][0] = matriz[i][2];
                matriz[i][columnas-1] = matriz[i][columnas-3];
            }
            //rellenar la primera y ultima fila para poder procesar los datos
            for(int i = 1; i<columnas; i++){
                matriz[0][i] = matriz[2][i];
                matriz[filas-1][i] = matriz[filas-3][i];
            }
            matriz[0][0] = matriz[2][2];
            matriz[0][columnas-1] = matriz[2][columnas-3];
            matriz[filas-1][0] = matriz[filas-3][2];
            matriz[filas-1][columnas-1] = matriz[filas-3][columnas-3];
            
            copiar_datos_bin(matriz, filas, columnas, fichero, "Sobel");

        }
    }else if(procesado == 2){//si se pide el histograma
        contador = histograma(cant_filas[myrank], columnas, matriz_aux);
        //se envian los datos del contador, se suman y guardan en el contador del root
        MPI_Reduce(contador, contador_aux, 256, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        if(myrank == 0){
            copiar_hist(contador_aux);
        }      
    }

    //liberar los espacios de memoria reservados
    for(int i = 0; i<cant_filas[myrank]+2; i++){
        free(matriz_aux[i]);
    }
    free(matriz_aux);
    free(cant_filas);
    free(pos_fila);

    

    if(myrank==0){
        //fin de medición de tiempo
        fin = MPI_Wtime();
        tiempo = (fin - inicio) * 1000;
        
        printf("\nTiempo medido: %.0f milisegundos.\n", tiempo);

        //liberar espacio de la matriz
        for(int i = 0; i<filas; i++){
            free(matriz[i]);
        }
        free(matriz);
    }

    MPI_Finalize();
    return;
}