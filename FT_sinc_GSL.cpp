#include <iostream>
#include <cmath>
#include <gsl/gsl_fft_complex.h>
#include <fstream>
#include "gnuplot.h"

#define REAL(z,i) (z[2*i])
#define IMAG(z,i) (z[2*i+1])

using namespace std;

double sinc(double x)
{
    if(x!=0)
        return sin(x)/x;
    else return 1.0;
}

double box(double k)
{
    if(k<=1 && k>=-1)
        return sqrt(M_PI/2);
    else return 0;
}

int main() {
    
    double pi = M_PI;
    
    int n=64;
    double x_min = -30.0;
    double x_max = 30.0;
    double dx = (x_max-x_min)/(n-1);
    double x_arr[n];
    double k_arr[n];
    double f_x[2*n];
    double f_q[2*n];
    
    int i;
    for(i=0;i<n;i++)
    {
        REAL(f_x,i) = sinc(x_min + i*dx);
        IMAG(f_x,i) = 0.0;
        x_arr[i] = x_min + i*dx;
        
        if(i<n/2)
            k_arr[i] = 2*pi*i/(dx*n);
        else k_arr[i] = 2*pi*i/(dx*n) - 2*pi/dx;
    }
    
    gsl_fft_complex_radix2_forward(f_x,1,n);
    
    ofstream data;
    data.open("a3.dat");
    
    for(i=0;i<n;i++)
    {
        REAL(f_q,i) = (dx/sqrt(2*pi))*(REAL(f_x,i)*cos(k_arr[i]*x_min) + IMAG(f_x,i)*sin(k_arr[i]*x_min));
        IMAG(f_q,i) = (dx/sqrt(2*pi))*(IMAG(f_x,i)*cos(k_arr[i]*x_min) - REAL(f_x,i)*sin(k_arr[i]*x_min));
        
        data << k_arr[i] << " " << REAL(f_q,i) << endl;
    }
    
    for(i=0;i<n;i++)
    {
        cout << REAL(f_q,i) << "   " << IMAG(f_q,i) << endl;
    }
    
    gnuplot q;
    q("set term postscript eps");
    q("set output \"ft_sinc_gsl.eps\" ");
    q("plot\'./a3.dat\' with points pointtype 5, ");
    
    return 0;
}
