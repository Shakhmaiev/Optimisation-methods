#include <iostream>
#include <iomanip>
#include <cmath>
using namespace std;

double F(double x1, double x2)
{
    return 4*x1*x1 + 3*x2*x2 + 2*x1*x2 - 18*x1 + 12*x2;
}

void grad(double x1, double x2, double &g1, double &g2)
{
    g1 = 8*x1 + 2*x2 - 18;
    g2 = 6*x2 + 2*x1 + 12;
}

int znak(double x)
{
    if (x > 0) return 1;
    if (x < 0) return -1;
    return 0;
}

void metod(int m, double h0)
{
    double x1 = 9;
    double x2 = 0;
    double eps = 0.01;
    double h = h0;

    double fx = F(x1, x2);

    int nf = 1;
    int ng = 0;
    int k = 0;

    if (m == 1)
        cout << "\nSferychna metryka\n\n";
    else if (m == 2)
        cout << "\nKubichna metryka\n\n";
    else
        cout << "\nOktaedrychna metryka\n\n";

    cout << left
         << setw(4)  << "k"
         << setw(11) << "x1"
         << setw(11) << "x2"
         << setw(13) << "F"
         << setw(11) << "g1"
         << setw(11) << "g2"
         << setw(12) << "|g|"
         << setw(11) << "d1"
         << setw(11) << "d2"
         << setw(10) << "h"
         << setw(8)  << "drob"
         << "Fnew"
         << endl;

    cout << string(125, '-') << endl;

    while (k < 5000)
    {
        double g1, g2;
        grad(x1, x2, g1, g2);
        ng++;

        double gn = sqrt(g1*g1 + g2*g2);

        if (gn <= eps)
            break;

        double d1, d2;

        // сферична
        if (m == 1)
        {
            d1 = -g1 / gn;
            d2 = -g2 / gn;
        }

        // кубічна
        else if (m == 2)
        {
            d1 = -znak(g1);
            d2 = -znak(g2);
        }

        // октаедрична
        else
        {
            d1 = 0;
            d2 = 0;

            if (abs(g1) >= abs(g2))
                d1 = -znak(g1);
            else
                d2 = -znak(g2);
        }

        int dr = 0;
        double xn1, xn2, fn;

        while (true)
        {
            xn1 = x1 + h*d1;
            xn2 = x2 + h*d2;

            fn = F(xn1, xn2);
            nf++;

            if (fn < fx)
                break;

            h = h / 2;
            dr++;

            if (h < 1e-10)
            {
                cout << "\nKrok zanadto malyy\n";
                return;
            }
        }

        cout << fixed << setprecision(6);

        cout << left
             << setw(4)  << k + 1
             << setw(11) << x1
             << setw(11) << x2
             << setw(13) << fx
             << setw(11) << g1
             << setw(11) << g2
             << setw(12) << gn
             << setw(11) << d1
             << setw(11) << d2
             << setw(10) << h
             << setw(8)  << dr
             << fn
             << endl;

        x1 = xn1;
        x2 = xn2;
        fx = fn;

        k++;
    }

    cout << "\nResult:\n";
    cout << "x* = (" << x1 << "; " << x2 << ")" << endl;
    cout << "F(x*) = " << fx << endl;
    cout << "Iteraciy: " << k << endl;
    cout << "NF = " << nf << endl;
    cout << "NgradF = " << ng << endl;
    cout << "Razom = " << nf + ng << endl;
}

int main()
{
    double h0 = 0.7;

    metod(1, h0);
    metod(2, h0);
    metod(3, h0);

    return 0;
}