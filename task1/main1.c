#include <stdio.h>
#include <math.h>
#include <stdlib.h>

int n[3] = {32, 64, 128};

double pi = 3.141592653589793238;

double u(double x, int i, int j){
    if(j==0)
        return sin(pi*x / 2);
    else if(j==1){
        return sin(pi*(n[i]-1)*x /2);
    }
    else
        return x*(1-x)*(1-x);
}

double func(double x, int i, int j){
    if(j == 0)
        return 2*x*x*sin(pi * x / 2) + (pi*pi / 4) * (1. + x) * sin(pi * x / 2) - (pi / 2)*cos(pi * x / 2);
    else if(j == 1){
        int N = n[i];
        return 2*x*x*sin(pi * (N-1) * x / 2) + (pi*pi*(N-1)*(N-1) / 4) * (1. + x)*sin(pi*(N-1)*x / 2) - (pi*(N-1) / 2)*cos(pi*(N-1)*x / 2);
    }
    else
        return 2*x*x*x*x*x - 4*x*x*x*x + 2*x*x*x - 9*x*x + 2*x + 3.0;
}

double *for_vec_b(int i, int j, double *vec_b){
    int k;
    double h = 1. / n[i];

    vec_b[0] = 0.;
    vec_b[1] = (h / 6) * (4*func(h, i, j) + func(2*h, i, j));
    vec_b[n[i]] = (h / 6) * (func((n[i]-1)*h, i, j) + 2*func(n[i]*h, i, j));

    for(k = 2; k < n[i]; k++)
        vec_b[k] = (h / 6) * (func((k-1)*h, i, j) + 4*func(k*h, i, j) + func((k+1)*h, i, j));

    return vec_b;
}

double *for_vec_b_new(int i, int j, double *vec_b){
    int k;
    double h = 1. / n[i];

    vec_b[0] = 0.;
    vec_b[1] = (h / 6) * (4*func(h, i, j) + func(2*h, i, j));
    vec_b[n[i]] = (h / 6) * (func((n[i]-1)*h, i, j) + 2*func(n[i]*h, i, j));

    for(k = 2; k < n[i]; k++)
        vec_b[k] = h * func((k+1)*h, i, j);

    return vec_b;
}

double *for_a(int i, double *a){
    int k;
    double h = 1. / n[i];

    a[0] = 0.;
    a[1] = 0.;

    for(k = 2; k < n[i]+1; k++)
        a[k] = (- h*h*h / 3)*((k-1)*(k-1) + k - 1 + 3./10) + k - 1./2 + 1./h;

    return a;
}

double *for_b(int i, double *b){
    int k;
    double h = 1. / n[i];

    b[0] = 0.;
    b[n[i]] = 0.;

    for(k = 1; k < n[i]; k++)
        b[k] = (- h*h*h / 3)*(k*k + k + 3./10) + k + 1./2 + 1./h;

    return b;
}

double *for_c(int i, double *c){
    int k;
    double h = 1. / n[i];

    c[0] = 0.;
    c[n[i]] = 2./h - 1./2 + (h/15)*(h*h - 5*h + 10);
    for(k = 1; k < n[i]; k++)
        c[k] = 2*k + 2./h + (h*h*h/15)*(20*k*k + 2);

    return c;

}

int main(){
    int i, j, k;
    double *vec_b, *a, *b, *c;
    double *ksi, *nu, *y, *x;
    double delta;

    for(j = 0; j < 3; j++) // итерация по функциям
        for(i = 0; i < 3; i++){     //итерация по N

            a = (double*)malloc((n[i]+1) * sizeof(double));
            b = (double*)malloc((n[i]+1) * sizeof(double));
            c = (double*)malloc((n[i]+1) * sizeof(double));

            a = for_a(i, a);
            b = for_b(i, b);
            c = for_c(i, c);

            vec_b = (double*)malloc((n[i]+1) * sizeof(double));

            if(j == 2)
                vec_b = for_vec_b_new(i, j, vec_b);
            else
                vec_b = for_vec_b(i, j, vec_b);

            ksi = (double*)malloc((n[i]+1) * sizeof(double));
            nu = (double*)malloc((n[i]+1) * sizeof(double));

            y = (double*)malloc((n[i]+1) * sizeof(double));
            x = (double*)malloc((n[i]+1) * sizeof(double));


            ksi[n[i]] = a[n[i]] / c[n[i]];
            for (k = n[i] - 1; k > 1 ; k--)
                ksi[k] = a[k] / (c[k] - b[k]*ksi[k+1]);


            nu[n[i]] = vec_b[n[i]] / c[n[i]];
            for (k = n[i] - 1; k > 0 ; k--)
                nu[k] = (vec_b[k] + b[k]*nu[k+1]) / (c[k] - b[k]*ksi[k+1]);

            y[0] = 0.;
            y[1] = nu[1];
            for(k = 1; k < n[i]; k++)
                y[k+1] = ksi[k+1]*y[k] + nu[k+1];

            x[0] = 0.;
            for(k = 1; k < n[i]+1; k++)
                x[k] = u((k* (1./ n[i])), i, j);

            delta = 0.0;

            for(k = 1; k < n[i]; k++)
                delta += ((x[k]-y[k])*(x[k]-y[k]));

            delta += ((x[n[i]]-y[n[i]])*(x[n[i]]-y[n[i]])) / 2;

            delta = delta / n[i];

            delta = sqrt(delta);

            printf("N = %d, nomer_func = %d, delta = %0.3e\n", n[i], j, delta);
        }

    return 0;
}