#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>

int numerico(char **argv){
    for(int j=2; j<5; j++){
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

void histograma(int filas, int columnas, unsigned char **matriz, FILE *f){

    //contar la cantidad de veces que se repite cada valor en la matriz
    int contador[256] = {0};
    for(int i = 1; i<filas+1; i++){
        for(int j = 1; j<columnas+1; j++){
            contador[matriz[i][j]]++;
        }
    }

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

void calculoMedia(unsigned char **matriz, int filas, int columnas, char *fichero){
    FILE *f;
    int suma;
    //asignar espacio para la nueva matriz de calculo de media
    unsigned char **matrizMedia;
    matrizMedia = (unsigned char**)malloc((filas) * sizeof(unsigned char*));
    for(int i = 0; i<filas; i++){
        matrizMedia[i] = (unsigned char*)malloc((columnas) * sizeof(unsigned char));
    }
    //calcular la media de cada dato de la matriz original y guardarlos en la matriz de medias
    for(int i = 1; i<filas+1; i++){
        for(int j = 1; j<columnas+1; j++){
            suma = 0;
            for(int k = 0; k<3; k++){
                for(int l = 0; l<3; l++){
                    //realizar la suma de los 9 elementos de la minimatriz
                    suma += matriz[i-1+k][j-1+l];
                }
            }
            //dividir la suma de la minimatriz entre 9 para hayar la media
            matrizMedia[i-1][j-1] = suma / 9;
        }
    }

    //copiar los datos en el fichero binario
    char newFile[25] = "Media";
    strcat(newFile,fichero);
    f = fopen(newFile, "wb");
    for(int i = 0; i<filas; i++){
        fwrite(matrizMedia[i], sizeof(unsigned char), columnas, f);
    }
    fclose(f);

    
    for(int i = 0; i<filas; i++){
        free(matrizMedia[i]);
    }
    free(matrizMedia);
    return;
}

void sobel(unsigned char **matriz, int filas, int columnas, char *fichero){
    FILE *fp;
    //asignar espacio para la matriz de resultados de C y F, en este caso son ints porque hay resultados que pueden ser negativos
    int c;
    int f;
    //declarar las matrices C y F
    int C[3][3] = {{-1, 0, 1},{-2, 0, 2},{-1, 0, 1}};
    int F[3][3] = {{-1, -2, -1},{0, 0, 0},{1, 2, 1}};
    //asignar espacio para la nueva matriz de calculo de media
    unsigned char **matrizSobel;
    matrizSobel = (unsigned char**)malloc((filas) * sizeof(unsigned char*));
    for(int i = 0; i<filas; i++){
        matrizSobel[i] = (unsigned char*)malloc((columnas) * sizeof(unsigned char));
    }
    //calcular la detección de bordes (sobel)  de cada dato de la matriz original y guardarlos en la matriz de medias
    for(int i = 2; i<filas; i++){
        for(int j = 2; j<columnas; j++){
            c = 0;
            f = 0;
            for(int k = 0; k<3; k++){
                for(int l = 0; l<3; l++){
                    c += matriz[i-1+k][j-1+l]*C[k][l];
                    f += matriz[i-1+k][j-1+l]*F[k][l];
                }
            }
            matrizSobel[i-1][j-1] = sqrt( c*c + f*f );
        }
    }
    //Extensión simétrica para tener una matriz resultante 512x512
    for(int i = 1; i<columnas; i++){
        matrizSobel[0][i] = matrizSobel[2][i];
        matrizSobel[filas-1][i] = matrizSobel[filas-3][i];
    }
    for(int i = 1; i<filas; i++){
        matrizSobel[i][0] = matrizSobel[i][2];
        matrizSobel[i][columnas-1] = matrizSobel[i][columnas-3];
    }
    matrizSobel[0][0] = matrizSobel[2][2];
	matrizSobel[0][columnas-1] = matrizSobel[2][columnas-3];
    matrizSobel[filas-1][0] = matrizSobel[filas-3][2];
	matrizSobel[filas-1][columnas-1] = matrizSobel[filas-3][columnas-3];

    //copiar los datos en el fichero binario
    char newFile[25] = "Sobel";
    strcat(newFile,fichero);
    fp = fopen(newFile, "wb");
    for(int i = 0; i<filas; i++){
        fwrite(matrizSobel[i], sizeof(unsigned char), columnas, fp);
    }
    fclose(fp);

    
    for(int i = 0; i<filas; i++){
        free(matrizSobel[i]);
    }
    free(matrizSobel);
    return;
}

void main(int argc, char *argv[]) {

    //comenzar a medir el tiempo que se demora
    clock_t inicio = clock();

    if(argc != 5){
        printf("Es necesario poner como argumento:\n 1- nombre del fichero(cadena de caracteres)\n 2- largo(numero)\n 3- ancho(numero) \n 4-la opción a escoger(numero)\n");
    }

    FILE *f;
    char *fichero = argv[1];
    if(numerico(argv) == 1){
        printf("Alguno de los valores que deberia ser numerico no lo es");
        return;
    }
    int filas = atoi(argv[2]);
    int columnas = atoi(argv[3]);
    unsigned char procesado = atoi(argv[4]);
    if(procesado < 1 || procesado > 3){
        printf("El argumento donde se escoge la opcion de procesado debe ser:\n1-Calculo de media\n2-Histograma\n3-Sobel");
        return;
    }
    unsigned char **matriz;
    
    //Abrimos el archivo binario .raw
    f = fopen(fichero, "rb");
    if (f == NULL) {
        printf("No se puede abrir el archivo binario\n");
        return;
    }

    //Le asignamos espacio de manera dinámica a la matriz en la que guardaremos los datos
    matriz = (unsigned char**)malloc((filas+2) * sizeof(unsigned char*));
    matriz[0] = (unsigned char*)malloc((columnas+2) * sizeof(unsigned char));
    for(int i = 1; i<=filas; i++){
        matriz[i] = (unsigned char*)malloc((columnas+2) * sizeof(unsigned char));
        //copiamos los datos del fichero en la matriz
        fread(&matriz[i][1], sizeof(unsigned char), columnas, f);
    }
    matriz[filas+1] = (unsigned char*)malloc((columnas+2) * sizeof(unsigned char));
    //cerramos el archivo binario
    fclose(f);


    //LLamar a función para hacer el histograma
    if(procesado == 2){
        histograma(filas, columnas, matriz, f);
    }

    //rellenar la primera y última fila para poder procesar los datos
    for(int i = 1; i<columnas+1; i++){
        matriz[0][i] = matriz[2][i];
        matriz[filas+1][i] = matriz[filas-1][i];
    }
    //Extención simétrica: rellenar la primera y última columna para poder procesar los datos
    for(int i = 1; i<filas+1; i++){
        matriz[i][0] = matriz[i][2];
        matriz[i][columnas+1] = matriz[i][columnas-1];
    }
    matriz[0][0] = matriz[2][2];
	matriz[0][columnas+1] = matriz[2][columnas-1];
    matriz[filas+1][0] = matriz[filas-1][2];
	matriz[filas+1][columnas+1] = matriz[filas-1][columnas-1];

    if(procesado == 1){
        calculoMedia(matriz, filas, columnas, fichero);
    }else if(procesado == 3){
        sobel(matriz, filas, columnas, fichero);
    }

    
    for(int i = 0; i<filas+2; i++){
        free(matriz[i]);
    }
    free(matriz);

    //fin de medición de tiempo
    clock_t fin = clock();
    
    double tiempo = (double)(fin - inicio) / CLOCKS_PER_SEC * 1000;
    
    printf("\nTiempo medido: %.0f milisegundos.\n", tiempo);


    return;
}