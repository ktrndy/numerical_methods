#include <stdio.h>
#include <math.h>
#include <stdlib.h>

int iter[3] = {10, 100, 1000};

double pi = 3.141592653589793238;

double u(double x, double t, int j){    // вроде ок
    if(j==0)
        return sqrt(30.)*t*x*(1. - x);
    else if(j==1)
        return sqrt(105.)*exp(-t + 1.)*x*x*(1. - x);
    else if(j == 2)
        return sqrt(630.)*t*x*x*(1. - x)*(1. - x);
    else
        return sqrt(2.)*exp(-pi*pi*(t - 1.))*sin(pi*x);
}

double phi_m(double x, double t, int j, int M){  // вроде ок
    if(j == 0)
        return sqrt(30.)*(x - x*x + 2*t) + (1. / (12*M*M))*(-2*sqrt(30.));
    else if(j == 1)
        return sqrt(105.)*exp(-t + 1.)*(x*x*x - x*x + 6*x - 2.) + (1. / (12*M*M))*sqrt(105.)*exp(-t + 1.)*(6*x - 2.);
    else if (j == 2)
        return sqrt(630.)*(x*x - 2*x*x*x + x*x*x*x - 2*t + 12*t*x - 12*t*x*x) + (1. / (12*M*M))*sqrt(630.)*(2. - 12*x + 12*x*x - 24*t);
    else
        return 0.0;
}

double *for_a(int N, int M, double *a){  //ок
    int k;

    a[0] = 0.;
    a[1] = 0.;

    for(k = 2; k < M; k++)
        a[k] = (1.*M*M)/2 - (1.*N)/12;

    return a;
}

double *for_b(int N, int M, double *b){  //ок
    int k;

    b[0] = 0.;
    b[M-1] = 0.;

    for(k = 1; k < M-1; k++)
        b[k] = (1.*M*M)/2 - (1.*N)/12;

    return b;
}

double *for_c(int N, int M, double *c){  //ок
    int k;

    c[0] = 0.;

    for(k = 1; k < M; k++)
        c[k] = (5. * N)/6 + 1.*M*M;

    return c;
}

double *for_f(int N, int M, double *f, int j, int n, const double *u_h){ //вроде ок
    int k;

    f[0] = 0.;

    f[1] = ((5.*N)/6 - 1.*M*M) * u_h[1] + ((1.*N)/12 + (1.*M*M)/2)*u_h[2] + phi_m(1./M, (1.*n + 1./2)/N, j, M);

    for(k = 2; k < M - 1; k++)
        f[k] = ((1.*N)/12 + (1.*M*M)/2)*u_h[k+1] + ((5.*N) / 6 - (1.*M*M)) * u_h[k] + ((1.*N)/12 + (1.*M*M)/2)*u_h[k-1] + phi_m(1.*k/M, (1.*n + 1./2)/N, j, M);

    f[M-1] = ((5.*N) / 6 - 1.*M*M) * u_h[M-1] + ((1.*N)/12 + (1.*M*M)/2)*u_h[M-2] + phi_m(1.*(M-1)/M, (1.*n + 1./2)/N, j, M);

    return f;
}

double *for_alpha(double *alpha, const double *a, const double *b, const double *c, int M){ //вроде ок
    int k;

    alpha[0] = 0.;
    alpha[1] = b[1] / c[1];

    for (k = 2; k < M-1; k++)
        alpha[k] = b[k] / (c[k] - a[k]*alpha[k-1]);


    return alpha;
}

double *for_beta(double *beta, const double *a, const double *c, const double *f, const double *alpha, int M){
    int k;

    beta[0] = 0.;
    beta[1] = f[1] / c[1];

    for (k = 2; k < M-1; k++)
        beta[k] = (f[k] + a[k] * beta[k-1]) / (c[k] - a[k]*alpha[k-1]);

    return beta;
}


int main(){
    int j, k, l, m, n;
    double *a, *b, *c, *f, *u_h;
    double *alpha, *beta, *u_1;
    double delta;

    for(j = 0; j < 4; j++) // итерация по функциям
            for(k = 0; k < 3; k++)  //по h  Mh = 1
                for (l = 0; l < 3; l++){ // по тау  Ntau = 1

                    int M = iter[k];
                    int N = iter[l];

                    a = (double*)malloc(M * sizeof(double));
                    b = (double*)malloc(M * sizeof(double));
                    c = (double*)malloc(M * sizeof(double));
                    f = (double*)malloc(M * sizeof(double));    // правая часть системы уравнений для прогонки

                    alpha = (double*)malloc((M-1) * sizeof(double));
                    beta = (double*)malloc((M-1) * sizeof(double));

                    a = for_a(N, M, a);
                    b = for_b(N, M, b);
                    c = for_c(N, M, c);

                    u_1 = (double*)malloc((M+1) * sizeof(double));
                    u_h = (double*)malloc((M+1) * sizeof(double));  // определяем начальный вектор, с него будем стартовать

                    u_h[0] = 0.;    // краевые для u^0 то есть при t = 0
                    u_h[M] = 0.;

                    for(m = 1; m < M; m++)
                        u_h[m] = u((1.*m) / M, 0.0, j);

                    alpha = for_alpha(alpha, a, b, c, M);

                    for(n = 0; n < N; n++){
                        f = for_f(N, M, f, j, n, u_h);
                        beta = for_beta(beta, a, c, f, alpha, M);

                        u_h[M-1] = (f[M-1] + a[M-1]*beta[M-2]) / (c[M-1] - a[M-1]*alpha[M-2]);
                        for(m = M-2; m > 0; m--)
                            u_h[m] = alpha[m]*u_h[m+1] + beta[m];   // вычислили u^(i+1)

                    }

                    u_1[0] = 0.;
                    u_1[M] = 0.;

                    for(m = 1; m < M; m++)
                        u_1[m] = u((1.*m)/M, 1., j);

                    delta = 0.;

                    for(m = 1; m < M; m++)
                        delta += (u_1[m]-u_h[m])*(u_1[m]-u_h[m]);

                    delta = delta / M;

                    delta = sqrt(delta);

                    printf("N = %d, M = %d, nomer_func = %d, delta = %0.3e\n", N, M, j, delta);
                }

    return 0;
}