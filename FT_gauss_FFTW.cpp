#include <iostream>
#include <cmath>
#include <fftw3.h>
#include <fstream>
#include "gnuplot.h"

using namespace std;

double gauss(double x)
{
    return exp(-x*x);
}

int main() {
    
    double pi = M_PI;
    
    int n=64;
    double x_min = -20.0;
    double x_max = 20.0;
    double dx = (x_max-x_min)/(n-1);
    double x_arr[n];
    double k_arr[n];
    fftw_complex w_p[n], tw_q[n], f_q[n];
    fftw_plan p;
    
    int i;
    for(i=0;i<n;i++)
    {
        w_p[i][0] = gauss(x_min + i*dx);
        w_p[i][1] = 0.0;
        x_arr[i] = x_min + i*dx;
        
        if(i<n/2)
            k_arr[i] = 2*pi*i/(dx*n);
        else k_arr[i] = 2*pi*i/(dx*n) - 2*pi/dx;
    }
    
    p = fftw_plan_dft_1d(n,w_p, tw_q, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(p);
    
    ofstream data;
    data.open("a2.dat");
    
    for(i=0;i<n;i++)
    {
        f_q[i][0] = (dx/sqrt(2*pi))*(tw_q[i][0]*cos(k_arr[i]*x_min) + tw_q[i][1]*sin(k_arr[i]*x_min));
        f_q[i][1] = (dx/sqrt(2*pi))*(tw_q[i][1]*cos(k_arr[i]*x_min) - tw_q[i][0]*sin(k_arr[i]*x_min));
        
        data << k_arr[i] << "  " << f_q[i][0] << endl;
    }
    data.close();
    
    for(i=0;i<n;i++)
    {
        cout << f_q[i][0] << "  " << f_q[i][1] << endl;
    }
    fftw_destroy_plan(p);
    
    gnuplot q;
    q("set term postscript eps");
    q("set output \"ft_gauss.eps\" ");
    q("plot\'./a2.dat\' u 1:2 w l");
    
    return 0;
}
