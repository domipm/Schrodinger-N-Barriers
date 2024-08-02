//PROGRAMA QUE CALCULA EL COEFICIENTE DE TRANSMISIÓN Y LA PROBABILIDAD DE ENCONTRAR LA PARTÍCULA A LA DERECHA

#include<cmath>
#include<iostream>
#include<iomanip>
#include"complex.h"
#include"gsl_rng.h"

#define MAX 10000 //Tamaño máximo matrices

#define PI 3.14159265358979323846 //Número PI
#define SEED 12824612670513 //Semilla números aleatorios

#define MAX_ITER 10000 //Número máximo de iteraciones a realizar en caso de no encontrar máximo

int N = 1000; //Subdivisiones espaciales eje X
int n_ciclos = 50;
int n_experimentos = 1000; //Número de experimentos a realizar

float k = 2*PI*n_ciclos/N;
float s = 1.0/(4.0*k*k);

float pot_centro = N/2.0; //Centro de la barrera
float pot_anchura = N/5.0; //Anchura de la barrera
float pot_altura = 0.5; //Parámetro lambda

float fonda_centro = N/4.0; //Posición media inicial de la partícula
float fonda_anchura = N/16.0; //Ancho de la función de onda inicial

float V[MAX]; //Vector potencial

fcomplex fonda[MAX]; //Función de onda
fcomplex A[MAX], alpha[MAX], B[MAX], beta[MAX], chi[MAX]; //Coeficientes
gsl_rng *tau;

using namespace std;

//Generamos potencial cuadrado centrado en (c) de anchura (a) y altura (h)
void sqr_pot(float c, float a, float h) {

    //Calculamos el potencial
    for (int j = 0; j < N; j++) {
        if ( (j > int(c-a/2.0) ) && (j < int(c+a/2.0)) ) V[j] = h*k*k;
        else V[j] = 0.0;
    }
    return;

}

//Probabilidad de detectar la partícula a la derecha
float pd() {

    float sum = 0.0;
    for (int j = int(pot_centro + pot_anchura/2.0 + 1); j < N; j++) sum += Cabs(Cmul(fonda[j], Conjg(fonda[j])));
    return sum;

}

//Genera onda plana gaussiana normalizada centrada en (c) de anchura (a)
void gaussian_normalized(float c, float a) {

    float norm = 0.0;
    for (int j = 0; j < N; j++) {
        if (j == 0 || j == N) fonda[j] = Complex(0.0,0.0);
        else fonda[j] = Cgauss(k*j, exp( -(j-c)*(j-c)/(2.0*a*a) ));
        norm += Cabs(Cmul(fonda[j], Conjg(fonda[j])));
    }
    for (int j = 0; j < N; j++) fonda[j] = Cdiv(fonda[j], Complex(sqrt((norm)),0));
    return;

}

//Función principal
int main() {

    extern gsl_rng *tau;
    tau = gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(tau, SEED);

    //Generamos potencial cuadrado y función de onda gaussiana
    sqr_pot(pot_centro, pot_anchura, pot_altura);
    gaussian_normalized(fonda_centro, fonda_anchura);

    //Calculamos vectores 'A' y 'alpha'
    for (int j = 0; j < N; j++) A[j] = Complex( -2-V[j], 2/s ); //Vector A_j^0
    alpha[N-1] = Complex(0.0, 0.0); //Condición de contorno alpha_N-1 = 0
    for (int j = N-2; j > 0; j--) alpha[j-1] = Cdiv( Complex(-1.0, 0.0), Cadd(A[j], alpha[j] ) ); //Vector alpha_j

    int mT = 0; //Detecciones a la derecha
    float prob_m = 0; //Probabilidad P_D(n_D)

    //Bucle sobre los experimentos
    for (int z = 1; z <= n_experimentos; z++) {

        cout << "Experimento " << z << endl; //Escribimos en pantalla el número del experimento realizado

        //Generamos función de onda inicial
        gaussian_normalized(fonda_centro, fonda_anchura);

        //Buscamos valor (t = nD) correspondiente al primer máximo local de P_D(t) mientras vamos evolucionando nD pasos
        float max = 0.0; //Valor máximo
        int cont = 0; //Contador de iteraciones
        while (true) {

            float norma = 0.0;

            //Vectores 'B' y 'beta'
            for (int j = 0; j < N; j++) B[j] = Cdiv( Cmul( Complex(0.0, 4.0), fonda[j]), Complex(s, 0.0) ); //Vector B_j
            beta[N-1] = Complex(0.0, 0.0); //Condición de contorno beta_N-1 = 0
            for (int j = N-2; j > 0; j--) beta[j-1] = Cdiv( Csub(B[j], beta[j]), Cadd(A[j], alpha[j]) ); //Vector beta_j

            //Vector 'chi'
            for (int j = 0; j < N-1; j++) chi[j+1] = Cadd( Cmul(alpha[j], chi[j]), beta[j] );

            //Calculamos la nueva función de onda
            for (int j = 0; j < N; j++) {
                fonda[j] = Csub( chi[j], fonda[j] );
                norma += Cabs(Cmul(fonda[j], Conjg(fonda[j])));
            }

            if ( pd() > max) max = pd(); //Encontramos el primer máximo local
            else if (pd() < max || cont >= MAX_ITER) break; //En cuanto comienza a disminuir, salimos del bucle

            cont++;

        }

        prob_m += max; //Probabilidad en el punto máximo
       
        //Simulamos medición: si (p < P_D) detectamos la partícula
        if (gsl_rng_uniform(tau) < pd()) mT += 1;

    }   

    //Calculamos coeficiente de transmisión (media de detecciones)
    float K = mT/float(n_experimentos);
    //Calculamos probabilidad de encontrar la partícula a la derecha
    prob_m /= n_experimentos;;

    //Mostramos los resultados en pantalla
    cout << "K = " << setprecision(15) << K << endl;
    cout << "P = " << setprecision(15) << prob_m << endl;

    return 0;

}