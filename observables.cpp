//PROGRAMA QUE CALCULA LOS VALORES ESPERADOS DE LOS OBSERVABLES Y REALIZA (O NO) LA MEDIDA
//OUTPUT: "observables.txt" fichero con los valores esperados de los operadores en función del tiempo

#include<cmath>
#include<iostream>
#include"complex.h"
#include"gsl_rng.h"

#define MAX 1000

#define PI 3.14159265358979323846
#define H 6.62607004e-34
#define MAX_ITER 1000
#define SEED 345678

int N = 1000; //Subdivisiones del eje espacial
int n_ciclos = 50;
int n_temp = 1000; //Número de iteraciones temporales a realizar

float k = 2*PI*n_ciclos/N;
float s = 1.0/(4.0*k*k);

bool medida = true; //Si medida==true, realizamos medida
bool maximizando = true; //Si maximizando==true, aplicamos el método de maximizar P_D obteniendo el tiempo de medida n_D
float t_medida = 500; //Cuando realizamos la medida (0, n_temp)

float pot_centro = N/2.0;
float pot_anchura = N/5.0;
float pot_altura = 0.5;

float fonda_centro = N/4.0;
float fonda_anchura = N/5.0;

float V[MAX]; //VECTOR POTENCIAL

//Vamos a calcular el valor esperado de los siguientes observables:
//Valore esperados: Posición (x), Momento (p), Energía cinética (ec), Energía potencial (ep), Energía total (ep)
float x, p, ec, ep, et;
//Incertidumbres en posición dx y momento dp
float dx, dp;
//Primeras y segundas derivadas de la función de onda
fcomplex d1[MAX], d2[MAX];

fcomplex fonda[MAX]; //FUNCIÓN DE ONDA (Complejo en general)
fcomplex A[MAX], alpha[MAX], B[MAX], beta[MAX], chi[MAX]; //COEFICIENTES
gsl_rng *tau;

using namespace std;

//Generamos potencial cuadrado centrado en (c) de anchura (a) y altura (h)
void sqr_pot(float c, float a, float h) {

    for (int j = 0; j < N; j++) {
        if ( (j > int(c-a/2.0) ) && (j < int(c+a/2.0)) ) V[j] = h*k*k;
        else V[j] = 0.0;
    }
    return;

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

//Probabilidad de detectar la partícula a la derecha
float pd() {

    float sum = 0.0;
    for (int j = int(pot_centro + pot_anchura/2.0 + 1); j < N; j++) sum += Cabs(Cmul(fonda[j], Conjg(fonda[j])));
    return sum;

}

//Calculamos el punto (0,N) para el cual la función de onda es máxima
int max_fonda() {

    float max = 0.0;
    int pos_max = 0;
    for (int j = 0; j < N; j++)
    if (Cabs(fonda[j]) > max) {
        max = Cabs(fonda[j]);
        pos_max = j;
    } 
    return pos_max;

}

//Valor esperado del operador posición <x>
float xe () {
    float x = 0;
    for (int j = 0; j < N; j++) x += j*Cabs(fonda[j])*Cabs(fonda[j]);
    return x;
}
//Valor esperado posición cuadrado <x^2>
float x2e () {
    float x2 = 0;
    for (int j = 0; j < N; j++) x2 += j*j*Cabs(fonda[j])*Cabs(fonda[j]);
    return x2;
}
//Valor esperado momento <p>
float pe() {
    fcomplex p = Complex(0, 0);
    for (int j = 0; j < N; j++) p = Cadd(p, Cmul( Complex(0, -1) , Cmul( Conjg(fonda[j]), d1[j]) ) );
    return p.r;
}
//Valor esperado momento cuadrado <p^2>
float p2e() {
    fcomplex p2 = Complex(0, 0);
    for (int j = 0; j < N; j++) p2 = Cadd(p2, Cmul( Conjg(fonda[j]), Cmul( Complex(-1,0), d2[j] ) ) );
    return p2.r;
}
//Valor esperado energía cinética
float ece() {
    fcomplex k = Complex(0, 0);
    for (int j = 0; j < N; j++) k = Cadd(k, Cmul( Complex(-1, 0), Cmul( Conjg(fonda[j]), d2[j] ) ) );
    return k.r;
}
//Valor esperado energía potencial
float epe() {
    float v = 0;
    for (int j = 0; j < N; j++) v += V[j]*Cabs(fonda[j])*Cabs(fonda[j]);
    return v;
}
//Valor esperado energía total
float ete() {
    fcomplex e = Complex(0, 0);
    for (int j = 0; j < N; j++) e = Cadd( e , Cadd( Cmul( Conjg(fonda[j]) , Cmul( Complex(-1,0) , d2[j] ) ) , Cmul( Conjg(fonda[j]) , Cmul( fonda[j], Complex(V[j], 0) ) ) ) );
    return e.r;
}

int main() {

    extern gsl_rng *tau;
    tau = gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(tau, SEED);

    FILE *out;
    out = fopen("observables.txt", "w"); //Fichero salida valor esperado observables

    fprintf(out, "t\t<x>\t\t\t\t<p>\t\t\t\t<V>\t\t\t\t<K>\t\t\t\t<E>\t\t\t\tdx\t\t\t\tdp\n");

    //Generamos potencial cuadrado y función de onda gaussiana
    sqr_pot(pot_centro, pot_anchura, pot_altura);
    gaussian_normalized(fonda_centro, fonda_anchura);

    //Calculamos vectores 'A' y 'alpha'
    for (int j = 0; j < N; j++) A[j] = Complex( -2-V[j], 2/s ); //VECTOR A_j^0
    alpha[N-1] = Complex(0.0, 0.0); //CONDICION CONTORNO alpha_N-1 = 0
    for (int j = N-2; j > 0; j--) alpha[j-1] = Cdiv( Complex(-1.0, 0.0), Cadd(A[j], alpha[j] ) ); //VECTOR alpha_j

    float max = 0.0; //Valor máximo de la probabilidad

    int t = 0;

    //Si maximizamos o no
    if (maximizando == true) {

        //Bucle hasta t = n_D
        while (true) {

            float norm = 0; //Igualamos la norma a cero en cada iteración

            //Vectores 'B' y 'beta'
            for (int j = 0; j < N; j++) B[j] = Cdiv( Cmul( Complex(0.0, 4.0), fonda[j]), Complex(s, 0.0) ); //VECTOR B_j
            beta[N-1] = Complex(0.0, 0.0); //CONDICIÓN DE CONTORNO beta_N-1 = 0
            for (int j = N-2; j > 0; j--) beta[j-1] = Cdiv( Csub(B[j], beta[j]), Cadd(A[j], alpha[j]) ); //VECTOR beta_j

            //Vector 'chi'
            for (int j = 0; j < N-1; j++) chi[j+1] = Cadd( Cmul(alpha[j], chi[j]), beta[j] );

            //Calculamos la nueva función de onda y su norma
            for (int j = 0; j < N; j++) {
                fonda[j] = Csub( chi[j], fonda[j] );
                norm += Cabs(Cmul(fonda[j], Conjg(fonda[j])));
            }

            if (pd() > max) max = pd(); //Encontramos el máximo
            else if (pd() < max || t >= MAX_ITER) break;

            //Calculamos derivada primera
            d1[0] = Cadd(fonda[1], Cmul(Complex(-1.0, 0.0), fonda[0]));
            for (int j = 1; j < N-1; j++) d1[j] = Cmul( Cadd(fonda[j+1], Cmul( Complex(-1.0, 0.0), fonda[j-1] )), Complex(0.5, 0.0) );
            d1[N-1] = Cadd(fonda[N-2], Cmul(Complex(-1.0, 0.0), fonda[N-1]) );

            //Calculamos derivada segunda
            d2[0] = Cadd(fonda[1], Cmul(Complex(-1.0, 0.0), d1[0]));
            for (int j = 1; j < N-1; j++) d2[j] = Cmul( Cadd(d1[j+1], Cmul( Complex(-1.0, 0.0), d1[j-1] )), Complex(0.5, 0.0) );
            d2[N-1] = Cadd(d1[N-2], Cmul(Complex(-1.0, 0.0), d1[N-1]) );

            //Calculamos los valores esperados
            x = xe();
            dx = sqrt(x2e() - x*x);
            p = pe();
            ec = ece();
            ep = epe();
            et = ete();
            dp = sqrt(p2e() - p*p);

            //Escribimos en el fichero los valores obtenidos para cada iteración
            fprintf(out, "%d\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\n", t, x, p, ep, ec, et, dx, dp);
            t++;

        }

        //Realizamos la medición en t = N_D
        if (medida == true) {

            //Calculamos derivada primera
            d1[0] = Cadd(fonda[1], Cmul(Complex(-1.0, 0.0), fonda[0]));
            for (int j = 1; j < N; j++) d1[j] = Cmul( Cadd(fonda[j+1], Cmul( Complex(-1.0, 0.0), fonda[j-1] )), Complex(0.5, 0.0) );
            d1[N-1] = Cadd(fonda[N-2], Cmul(Complex(-1.0, 0.0), fonda[N-1]) );

            //Calculamos derivada segunda
            d2[0] = Cadd(fonda[1], Cmul(Complex(-1.0, 0.0), d1[0]));
            for (int j = 1; j < N-1; j++) d2[j] = Cmul( Cadd(d1[j+1], Cmul( Complex(-1.0, 0.0), d1[j-1] )), Complex(0.5, 0.0) );
            d2[N-1] = Cadd(d1[N-2], Cmul(Complex(-1.0, 0.0), d1[N-1]) );

            //Calculamos los valores esperados
            x = xe();
            dx = sqrt(x2e() - x*x);
            p = pe();
            ec = ece();
            ep = epe();
            et = ete();
            dp = sqrt(p2e() - p*p);

            //Escribimos en el fichero los valores obtenidos para cada iteración
            fprintf(out, "%d\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\n", t, x, p, ep, ec, et, dx, dp);

            float prob = 0.0;
            if (gsl_rng_uniform(tau) < pd()) { //Detectamos la partícula a la derecha
                for (int j = 0; j < pot_centro; j++) {
                    prob += Cabs(fonda[j])*Cabs(fonda[j]);
                    fonda[j] = Complex(0.0, 0.0);
                } 
                for (int j = pot_centro; j < N; j++) fonda[j] = Cdiv( fonda[j], Complex(1-prob, 0) );
            } else { //No detectamos la partícula a la derecha
                for (int j = pot_centro; j < N; j++) {
                    prob += Cabs(fonda[j])*Cabs(fonda[j]);
                    fonda[j] = Complex(0.0, 0.0); 
                    } 
                for (int j = 0; j < pot_centro - pot_anchura/2; j++) fonda[j] = Cdiv( fonda[j], Complex(1-prob, 0) );
            }

            //Escribimos en pantalla el valor esperado de cada observable al realizar la medición
            cout << "Medida en t = " << t << endl;
            cout << "<x> = " << x << endl;
            cout << "dx = " << dx << endl;
            cout << "<p> = " << p << endl;
            cout << "dp = " << dp << endl;
            cout << "<Ec> = " << ec << endl;
            cout << "<Ep> = " << ep << endl;
            cout << "<Et> = " << et << endl;

        }

        //Dejamos evolucionar el sistema los pasos restantes
        while (t < n_temp) {

            float norm = 0.0;

            //Vectores 'B' y 'beta'
            for (int j = 0; j < N; j++) B[j] = Cdiv( Cmul( Complex(0.0, 4.0), fonda[j]), Complex(s, 0.0) ); //VECTOR B_j
            beta[N-1] = Complex(0.0, 0.0); //CONDICIÓN DE CONTORNO beta_N-1 = 0
            for (int j = N-2; j > 0; j--) beta[j-1] = Cdiv( Csub(B[j], beta[j]), Cadd(A[j], alpha[j]) ); //VECTOR beta_j

            //Vector 'chi'
            for (int j = 0; j < N-1; j++) chi[j+1] = Cadd( Cmul(alpha[j], chi[j]), beta[j] );

            //Calculamos la nueva función de onda y su norma
            for (int j = 0; j < N; j++) {
                fonda[j] = Csub( chi[j], fonda[j] );
                norm += Cabs(Cmul(fonda[j], Conjg(fonda[j])));
            }

            //Calculamos derivada primera
            d1[0] = Cadd(fonda[1], Cmul(Complex(-1.0, 0.0), fonda[0]));
            for (int j = 1; j < N-1; j++) d1[j] = Cmul( Cadd(fonda[j+1], Cmul( Complex(-1.0, 0.0), fonda[j-1] )), Complex(0.5, 0.0) );
            d1[N-1] = Cadd(fonda[N-2], Cmul(Complex(-1.0, 0.0), fonda[N-1]) );

            //Calculamos derivada segunda
            d2[0] = Cadd(fonda[1], Cmul(Complex(-1.0, 0.0), d1[0]));
            for (int j = 1; j < N-1; j++) d2[j] = Cmul( Cadd(d1[j+1], Cmul( Complex(-1.0, 0.0), d1[j-1] )), Complex(0.5, 0.0) );
            d2[N-1] = Cadd(d1[N-2], Cmul(Complex(-1.0, 0.0), d1[N-1]) );

            //Calculamos los valores esperados
            x = xe();
            dx = sqrt(x2e() - x*x);
            p = pe();
            ec = ece();
            ep = epe();
            et = ete();
            dp = sqrt(p2e() - p*p);

            //Escribimos en el fichero los valores obtenidos para cada iteración
            fprintf(out, "%d\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\n", t, x, p, ep, ec, et, dx, dp);

            t++;

        }

    } else {

        //Bucle principal
        for (int t = 0; t < n_temp; t++) {

            float norm = 0; //IGUALAMOS LA NORMA A CERO EN CADA ITERACION

            //Vectores 'B' y 'beta'
            for (int j = 0; j < N; j++) B[j] = Cdiv( Cmul( Complex(0.0, 4.0), fonda[j]), Complex(s, 0.0) ); //VECTOR B_j
            beta[N-1] = Complex(0.0, 0.0); //CONDICIÓN DE CONTORNO beta_N-1 = 0
            for (int j = N-2; j > 0; j--) beta[j-1] = Cdiv( Csub(B[j], beta[j]), Cadd(A[j], alpha[j]) ); //VECTOR beta_j

            //Vector 'chi'
            for (int j = 0; j < N-1; j++) chi[j+1] = Cadd( Cmul(alpha[j], chi[j]), beta[j] );

            //Calculamos la nueva función de onda y su norma
            for (int j = 0; j < N; j++) {
                fonda[j] = Csub( chi[j], fonda[j] );
                norm += Cabs(Cmul(fonda[j], Conjg(fonda[j])));
            }

            //Calculamos derivada primera
            d1[0] = Csub(fonda[1], fonda[0]);
            for (int j = 1; j < N-1; j++) d1[j] = Cmul( Csub(fonda[j+1], fonda[j-1]), Complex(0.5, 0.0) );
            d1[N-1] = Csub(fonda[N-2], fonda[N-1]);

            //Calculamos derivada segunda
            d2[0] = Csub(fonda[1], d1[0]);
            for (int j = 1; j < N-1; j++) d2[j] = Cmul( Csub(d1[j+1], d1[j-1]), Complex(0.5, 0.0) );
            d2[N-1] = Csub(d1[N-2], d1[N-1]) ;

            //Calculamos los valores esperados
            x = xe();
            dx = sqrt(x2e() - x*x);
            p = pe();
            ec = ece();
            ep = epe();
            et = ete();
            dp = sqrt(p2e() - p*p);

            //Escribimos en el fichero los valores obtenidos para cada iteración
            fprintf(out, "%d\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\n", t, x, p, ep, ec, et, dx, dp);

            //Realizamos una medición tras un tiempo
            if (t == t_medida && medida == true) {
                float prob = 0.0;
                if (gsl_rng_uniform(tau) < pd()) { //Detectamos la partícula a la derecha
                    for (int j = 0; j < pot_centro; j++) {
                        prob += Cabs(fonda[j])*Cabs(fonda[j]);
                        fonda[j] = Complex(0.0, 0.0);
                    } 
                    for (int j = pot_centro; j < N; j++) fonda[j] = Cdiv( fonda[j], Complex(1-prob, 0) );
                } else { //No detectamos la partícula a la derecha
                    for (int j = pot_centro; j < N; j++) {
                        prob += Cabs(fonda[j])*Cabs(fonda[j]);
                        fonda[j] = Complex(0.0, 0.0); 
                    } 
                    for (int j = 0; j < pot_centro - pot_anchura/2; j++) fonda[j] = Cdiv( fonda[j], Complex(1-prob, 0) );
                }

                //Escribimos en pantalla el valor esperado de cada observable al realizar la medición
                cout << "Medida en t = " << t << endl;
                cout << "<x> = " << x << endl;
                cout << "dx = " << dx << endl;
                cout << "<p> = " << p << endl;
                cout << "dp = " << dp << endl;
                cout << "<Ec> = " << ec << endl;
                cout << "<Ep> = " << ep << endl;
                cout << "<Et> = " << et << endl;

            }

        }

    }

    fclose(out); //Cerramos el fichero de salida

    return 0;

}