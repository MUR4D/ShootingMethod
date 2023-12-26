
#include <stdio.h>
#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include<iostream>
#include<string>
#define delta 1.e-6
#define N 2
#define T 1.0
using namespace std;


void fcn(double t, double z[N], double f[N])
{

    f[0] = z[1];
    f[1] = 2 * z[1] + 3 * z[0] + exp(4 * t);

}

double dmax4(double ay, double ay1, double ay2, double ay3)
{
    double w;
    w = ay;
    if (ay1 > w) w = ay1;
    if (ay2 > w) w = ay2;
    if (ay3 > w) w = ay3;
    return w;
}

double dmax1(double a, double b)
{
    if (a > b) { return a; }
    return b;
}
double dmin1(double a, double b)
{
    if (a < b) { return a; }
    return b;
}
//константы в методе Дормана-Принса, числа взяты из книги Хайрер, Нерсетт, Ваннер "Решение ОДУ"

constexpr double  c2 = 1. / 5.,
c3 = 3. / 10.,
c4 = 4. / 5.,
c5 = 8. / 9.,
c6 = 1.,
c7 = 1.,
a21 = 1. / 5.,
a31 = 3. / 40.,
a32 = 9. / 40.,
a41 = 44. / 45.,
a42 = -56. / 15.,
a43 = 32. / 9.,
a51 = 19372. / 6561.,
a52 = -25360. / 2187.,
a53 = 64448. / 6561.,
a54 = -212. / 729.,
a61 = 9017. / 3168.,
a62 = -355. / 33.,
a63 = 46732. / 5247.,
a64 = 49. / 176.,
a65 = -5103. / 18656.,
a71 = 35. / 384.,
a72 = 0,
a73 = 500. / 1113.,
a74 = 125. / 192.,
a75 = -2187. / 6784.,
a76 = 11. / 84.,
b1 = 35. / 384.,
b2 = 0,
b3 = 500. / 1113.,
b4 = 125. / 192.,
b5 = -2187. / 6784.,
b6 = 11. / 84.,
b7 = 0,
r1 = 5179. / 57600.,
r2 = 0,
r3 = 7571. / 16695.,
r4 = 393. / 640.,
r5 = -92097. / 339200.,
r6 = 187. / 2100.,
r7 = 1. / 40.;


double dopri5(double* z, double t0, double tend, double eps) {
    int i = 0, reject = 0, nacc = 0, nrej = 0;
    double k1[N]{}, k2[N]{}, k3[N]{}, k4[N]{}, k5[N]{}, k6[N]{}, k7[N]{}, z1[N]{}, b[N]{}, r[N]{};
    double err, uround = 2.e-16, h = 0.1, denon, fac, tph, lambda = 0.5, hnew, dglob = 0.0;
    double hmax = 1.57, t = 0;
   
    int nsteps = 0;
    int counter = 0;
    while (eps < (tend - t) && eps * fabs(tend) < (tend - t)) {
        if ((t + h - tend) > 0.e0) h = tend - t;
        fcn(t, z, k1);
        for (i = 0; i < (N); i++) z1[i] = z[i] + h * a21 * k1[i];
        fcn(t + h * c2, z1, k2);
        for (i = 0; i < (N); i++) z1[i] = z[i] + h * (a31 * k1[i] + a32 * k2[i]);
        fcn(t + h * c3, z1, k3);
        for (i = 0; i < (N); i++) z1[i] = z[i] + h * (a41 * k1[i] + a42 * k2[i] + a43 * k3[i]);
        fcn(t + h * c4, z1, k4);
        for (i = 0; i < (N); i++) z1[i] = z[i] + h * (a51 * k1[i] + a52 * k2[i] + a53 * k3[i] + a54 * k4[i]);
        fcn(t + h * c5, z1, k5);
        for (i = 0; i < (N); i++) z1[i] = z[i] + h * (a61 * k1[i] + a62 * k2[i] + a63 * k3[i] + a64 * k4[i] + a65 * k5[i]);
        fcn(t + h * c6, z1, k6);
        for (i = 0; i < (N); i++) z1[i] = z[i] + h * (a71 * k1[i] + a72 * k2[i] + a73 * k3[i] + a74 * k4[i] + a75 * k5[i] + a76 * k6[i]);
        fcn(t + h * c7, z1, k7);
        for (i = 0; i < (N); i++) {
            b[i] = z[i] + h * (b1 * k1[i] + b2 * k2[i] + b3 * k3[i] + b4 * k4[i] + b5 * k5[i] + b6 * k6[i] + b7 * k7[i]);
            r[i] = z[i] + h * (r1 * k1[i] + r2 * k2[i] + r3 * k3[i] + r4 * k4[i] + r5 * k5[i] + r6 * k6[i] + r7 * k7[i]);
        }
        tph = t + h;
        err = 0.0;
        for (i = 0; i < (N); i++) {
            denon = dmax4(1.e-6, fabs(b[i]), fabs(z[i]), 2.e0 * uround / eps);
            err += ((b[i] - r[i]) / denon) * ((b[i] - r[i]) / denon);
        }
        err = sqrt(err / N);
        fac = dmax1((1.e0 / 6.e0), dmin1(3.e0, pow(err / eps, 0.125) / .9e0));

        hnew = h / fac;
        if (err <= eps)
        {
            nacc++;
            nsteps++;
            for (i = 0; i < (N); i++) { z[i] = b[i]; }
            t = tph;
            if (fabs(hnew) > hmax) hnew = hmax;
            if (reject) hnew = dmin1(fabs(hnew), fabs(h));
            reject = 0;
            lambda = -0.5 + 1.5 * (z[0] * z[0]);
            dglob = err + dglob * exp(lambda * h);
            h = hnew;
        }
        else
        {
            nrej++;
            nsteps++;
            reject = 1;
            h = hnew;
        }
        counter++;
        std::cout  << "y(" << t << ") = " << z[0] << std::endl;
    }

    //printf("nacc:%d\n", nacc);
    //printf("nrej:%d\n", nrej);
    //if( < eps)
    return z[0];//abs(z[1] + z[0] - 4 * exp(3) - exp(4));
}
// scanf("%f%f%f%f%f", &t0, &x0, &tn, &xn, &h);
 /*system("cls");
 printf("%f", t0);
 system("pause");*/

double Solution(double x)
{

    return -exp(-x) + 3*exp(3 * x) + 0.8 * exp(4 * x);
}


int main(int argc, char* argv[]) {

    int i;
    double z[N], tend = 1.0, eps, x1, y1, glob1;
    eps = 1.e-5;
    double res;
    double t0 = 0, x0 = 0, h = std::stod(argv[3]), tn = 1.0, xn = 0, z0, m1, m2, m3, b, b1, b2, b3, e;
    double h1 = 0;
    //  printf("\n  Enter t0,x0,tn,xn:");
      //std::cin >> t0 >> x0 >> tn >> xn;
    /*printf("\n  Enter the trial M1:");
    std::cin >> m1;*/
  //  b = xn;
  
    
        z[1] = i;
        z[0] = i-0.6;
        res = dopri5(z, t0, T, eps);
        if (fabs(res - 31.3731 <  eps))
        {
            cout << "y (" << h1 << ") = " << res << endl;
        }h1 += h;
        //cout << res << endl;
        //printf("\ny is %f\n", b1);
    
   // cout << res
    
  //  printf("\nB1 is %f", Solution(0));
    
    return 0;
}


//glob1 = dopri5(z, t0, T, eps);
/*printf("%f\n", z[1]);*/
/*for (i = 0; i < N; i++)
{
    printf("%.8lf \t", z[i]);
}*/
