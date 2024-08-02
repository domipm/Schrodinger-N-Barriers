//PROGRAMA QUE CALCULA EL COEFICIENTE DE TRANSMISIÓN Y LA PROBABILIDAD DE ENCONTRAR
//LA PARTÍCULA A LA DERECHA PARA VARIAS BARRERAS DE POTENCIAL

#include<cmath>
#include<iostream>
#include<iomanip>
#include"gsl_rng.h"
#include"complex.h"

#define MAX 1000 //Tamaño máximo de vectores

#define PI 3.14159265358979323846
#define SEED 12824612670513

int N = 1000; //Discretización eje X
int n_ciclos = 50;
int n_experimentos = 1000; //Número de experimentos a realizar
int n_barreras = 4; //Número de barreras de potencial

float k = 2*PI*n_ciclos/N;
float s = 1.0/(4.0*k*k);

float pot_anchura = N/16.0; //Anchura de cada barrera de potencial
float pot_altura = 5; //Parámetro lambda
float pot_distancia = N/5.0; //Distancia entre cada barrera de potencial

float fonda_centro = 0.0; //Centro función de onda inicial
float fonda_anchura = N/64.0; //Anchura función de onda inicial

float V[MAX]; //Vector potencial

fcomplex fonda[MAX]; //Función de onda
fcomplex A[MAX], alpha[MAX], B[MAX], beta[MAX], chi[MAX]; //Coeficientes A, alpha, B, beta, chi
gsl_rng *tau;

using namespace std;

//Generamos (n) barreras de potencial de anchura (a) y altura (h) equiespaciadas una distancia (d)
void sqr_pot(int n, float a, float h, float d) {

    //Calculamos el centro de la primera barrera
    int c = N / (n+1);
    //Calculamos el potencial
    for (int j = 0; j < N; j++) V[j] = 0.0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < N; j++) {
            if ( (j > int( (c+c*i) - a/2.0 ) ) && (j < int( (c+c*i) + a/2.0 )) ) V[j] = h*k*k;
        }
    }
    return;

}

//Probabilidad de detectar la partícula a la derecha de todas las barreras
float pd() {

    //Calculamos extremo derecho de la última barrera de potencial
    int c = N / (n_barreras + 1) + N / (n_barreras + 1) * (n_barreras - 1) + 1 + pot_anchura / 2;
    //Calculamos la probabilidad
    float sum = 0.0;
    for (int j = c; j < N; j++) sum += Cabs(Cmul(fonda[j], Conjg(fonda[j])));
    return sum;

}

//Genera onda plana gaussiana normalizada centrada en (c) de anchura (a)
void gaussian_normalized(float c, float a) {

    float norm = 0.0;
    //Calculamos la función de onda y su norma
    for (int j = 0; j < N; j++) {
        if (j == 0 || j == N) fonda[j] = Complex(0.0,0.0);
        else fonda[j] = Cgauss(k*j, exp( -(j-c)*(j-c)/(2.0*a*a) ));
        norm += Cabs(Cmul(fonda[j], Conjg(fonda[j])));
    }
    //Normalizamos la función de onda
    for (int j = 0; j < N; j++) fonda[j] = Cdiv(fonda[j], Complex(sqrt((norm)),0));
    return;

}

int main() {

    extern gsl_rng *tau;
    tau = gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(tau, SEED);

    //Generamos función de onda inicial
    gaussian_normalized(fonda_centro, fonda_anchura);
    //Generamos potencial cuadrado y función de onda gaussiana
    sqr_pot(n_barreras, pot_anchura, pot_altura, pot_distancia);

    //Calculamos vectores 'A' y 'alpha'
    for (int j = 0; j < N; j++) A[j] = Complex( -2-V[j], 2/s ); //Vector A_j^0
    alpha[N-1] = Complex(0.0, 0.0); //Condición de contorno alpha_N-1 = 0
    for (int j = N-2; j > 0; j--) alpha[j-1] = Cdiv( Complex(-1.0, 0.0), Cadd(A[j], alpha[j] ) ); //Vector alpha_j

    int mT = 0; //Detecciones a la derecha
    float prob_m = 0; //Probabilidad P_D(n_D) promediada

    for (int z = 1; z <= n_experimentos; z++) {

        cout << "Experimento " << z << endl;

        //Generamos función de onda 
        gaussian_normalized(fonda_centro, fonda_anchura);

        //Buscamos valor (t=nD) correspondiente al primer máximo local de P_D(t)
        float max = 0.0;
        int cont = 0;
        //Buscamos valor (t = n_D) primer máximo local de P_D(t)
        while (true) {

            //Vectores 'B' y 'beta'
            for (int j = 0; j < N; j++) B[j] = Cdiv( Cmul( Complex(0.0, 4.0), fonda[j]), Complex(s, 0.0) ); //Vector B_j
            beta[N-1] = Complex(0.0, 0.0); //Condición de contorno beta_N-1 = 0
            for (int j = N-2; j > 0; j--) beta[j-1] = Cdiv( Csub(B[j], beta[j]), Cadd(A[j], alpha[j]) ); //Vector beta_j

            //Vector 'chi'
            for (int j = 0; j < N-1; j++) chi[j+1] = Cadd( Cmul(alpha[j], chi[j]), beta[j] );

            //Calculamos la nueva función de onda
            for (int j = 0; j < N; j++) {
                fonda[j] = Csub( chi[j], fonda[j] );
            }

            if ( pd() > max ) max = pd();
            else if (pd() < max) break;

            cont++;

        }

        prob_m += max;

        //Si (p > P_D) detectamos la partícula
        if (gsl_rng_uniform(tau) < pd()) mT += 1;

    }   

    //Coeficiente de transmisión
    float K = float(mT) / float(n_experimentos);
    //Probabilidad en el punto máximo
    prob_m /= float(n_experimentos);

    cout << "K = " << K << endl;
    cout << "P = " << prob_m << endl;

    fcloseall();

    return 0;

}