#include<cmath>
#include<iostream>
#include"complex.h"

#define MAX 1000

#define PI 3.14159265358979323846
#define H 6.62607004e-34

int N = 1000, n_ciclos = 50, n_temp;
double lam, k, s, norm = 0;
double V[MAX]; //VECTOR POTENCIAL

int n_barreras = 4; //Número de barreras de potencial
float pot_anchura = N/16.0; //ORIGINALMENTE N/5.0
float pot_altura = 1.25; //Parámetro lambda
float pot_distancia = N/5.0;

fcomplex fonda[MAX]; //FUNCIÓN DE ONDA (Complejo en general)
fcomplex A[MAX], alpha[MAX], B[MAX], beta[MAX], chi[MAX]; //COEFICIENTES
fcomplex i = Complex(0,1);

using namespace std;

int main() {


    // ? F I C H E R O  S A L I D A //

    FILE *out;
    out = fopen("fonda.dat", "w"); //FICHERO SALIDA MÓDULO FUNCIÓN DE ONDA
    FILE *pot;
    pot = fopen("pot.dat", "w"); //FICHERO SALIDA POTENCIAL
    FILE *norma;
    norma = fopen("norm.dat", "w"); //FICHERO SALIDA NORMAL DE LA FUNCIÓN DE ONDA


    // ? P A R A M E T R O S  I N I C I A L E S //

    N = 1000; //DIVISIONES DEL EJE X
    n_ciclos = 100; //CICLOS A REALIZAR
    n_temp = 10000; //NUMERO DE ITERACIONES TEMPORALES
    lam = 0.8; //PARAMETRO LAMBDA (Altura del potencial)

    k = 2*PI*double(n_ciclos)/double(N); //FACTOR K REESCALADO
    s = 1/(4*k*k); //ESPACIADO TEMPORAL REESCALADO


    // ? P O T E N C I A L //


        //Calculamos el centro de la primera barrera
    int c = N / (n_barreras+1);

    //Calculamos el potencial
    for (int j = 0; j < N; j++) V[j] = 0.0;

    for (int i = 0; i < n_barreras; i++) {
        for (int j = 0; j < N; j++) {
            if ( (j > int( (c+c*i) - pot_anchura/2.0 ) ) && (j < int( (c+c*i) + pot_anchura/2.0 )) ) V[j] = pot_altura*k*k;
        }
    }

    for (int j = 0; j < N; j++) fprintf(pot, "%.15lf\n", V[j]); //Escribimos potencial en fichero


    // ? F U N C I Ó N  D E  O N D A  I N I C I A L //


    // * C Á L C U L O  D E  L A  F U N C I Ó N //

    for (int j = 0; j < N; j++) {
        if (j == 0 || j == (N)) fonda[j] = Complex(0,0); //CONDICIONES DE CONTORNO
        else {
            fonda[j] = Complex( cos(k*j) * exp( ( -8*double(4*j-N)*double(4*j-N) ) / double(N*N) ) , sin(k*j) * exp( ( -8*double(4*j-N)*double(4*j-N) ) / double(N*N) ) ); //CALCULAMOS FUNCIÓN DE ONDA
            norm += Cabs(Cmul(fonda[j], Conjg(fonda[j]))); //CALCULAMOS LA NORMA DE LA FUNCIÓN DE ONDA
        }
    }

    // * N O R M A L I Z A C I O N //
    
    for (int j = 0; j < N; j++) {
        fonda[j] = Cdiv(fonda[j], Complex(sqrt((norm)),0)); //DIVIDIMOS ENTRE NORMA
        fprintf(out, "%.15lf\n", Cabs(fonda[j])); //ESCRIBIMOS LA FUNCIÓN DE ONDA YA NORMALIZADA EN CADA PUNTO j DEL EJE x
    }
    fprintf(out, "\n"); //LINEA NUEVA PARA SEPARAR LAS ITERACIONES
    

    // * V E C T O R E S  'A_j^0'  Y  'A L P H A' //

    for (int j = 0; j < N; j++) A[j] = Complex( -2-V[j], 2/s ); //VECTOR A_j^0
    alpha[N-1] = Complex(0.0, 0.0); //CONDICION CONTORNO alpha_N-1 = 0
    for (int j = N-2; j > 0; j--) alpha[j-1] = Cdiv( Complex(-1.0, 0.0), Cadd(A[j], alpha[j] ) ); //VECTOR alpha_j


    // ? A L G O R I T M O //


    for (int t = 0; t < n_temp; t++) {

        norm = 0; //IGUALAMOS LA NORMA A CERO EN CADA ITERACION

        // * V E C T O R E S  'B_j'  Y  'B E T A' //

        for (int j = 0; j < N; j++) B[j] = Cdiv( Cmul( Complex(0.0, 4.0), fonda[j]), Complex(s, 0.0) ); //VECTOR B_j
        beta[N-1] = Complex(0.0, 0.0); //CONDICIÓN DE CONTORNO beta_N-1 = 0
        for (int j = N-2; j > 0; j--) beta[j-1] = Cdiv( Csub(B[j], beta[j]), Cadd(A[j], alpha[j]) ); //VECTOR beta_j

        // * V E C T O R  'C H I' //

        for (int j = 0; j < N-1; j++) chi[j+1] = Cadd( Cmul(alpha[j], chi[j]), beta[j] );

        // * F U N C I Ó N  D E  O N D A //

        for (int j = 0; j < N; j++) {
            fonda[j] = Csub( chi[j], fonda[j] );
            norm += Cabs(Cmul(fonda[j], Conjg(fonda[j])));
        }

        fprintf(norma, "%d\t%.15lf\n", t, sqrt(norm)); //ESCRIBIMOS EN EL FICHERO LA RAIZ DE LA NORMA AL CUADRADO (EN CADA ITERACION)
        for (int j = 0; j < N; j++) fprintf(out, "%.15lf\n", Cabs(fonda[j])); //ESCRIBIMOS LA FUNCIÓN DE ONDA EN EL FICHERO (SU MÓDULO)
        fprintf(out, "\n"); //LINEA NUEVA TRAS CADA ITERACIÖN PARA SEPARAR LAS FUNCIONES DE ONDA

    }

    fclose(out);
    fclose(pot);
    fclose(norma);

    return 0;

}