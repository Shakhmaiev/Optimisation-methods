#include <iostream>
#include <iomanip>
#include <cmath>
using namespace std;

int nf = 0;
int ng = 0;

double F(double x1, double x2)
{
    nf++;

    return 4*x1*x1 + 3*x2*x2
         + 2*x1*x2 - 18*x1 + 12*x2;
}

void grad(double x1, double x2, double &g1, double &g2)
{
    ng++;

    g1 = 8*x1 + 2*x2 - 18;
    g2 = 2*x1 + 6*x2 + 12;
}

double phi(double h,
           double x1, double x2,
           double d1, double d2)
{
    return F(x1 + h*d1, x2 + h*d2);
}

void svenn(double x1, double x2,
           double d1, double d2,
           double &a, double &b)
{
    double del = 0.01;

    double h0 = 0;
    double f0 = phi(h0, x1, x2, d1, d2);

    double h1 = del;
    double f1 = phi(h1, x1, x2, d1, d2);

    if (f1 >= f0)
    {
        h1 = -del;
        f1 = phi(h1, x1, x2, d1, d2);

        if (f1 >= f0)
        {
            a = -del;
            b = del;
            return;
        }

        del = -del;
    }

    double hp = h0;
    double hc = h1;

    int k = 1;

    while (true)
    {
        double hn = hc + pow(2.0, k) * del;

        double fc = phi(hc, x1, x2, d1, d2);
        double fn = phi(hn, x1, x2, d1, d2);

        if (fn >= fc)
        {
            a = min(hp, hn);
            b = max(hp, hn);
            return;
        }

        hp = hc;
        hc = hn;
        k++;
    }
}

double golden(double a, double b,
              double x1, double x2,
              double d1, double d2)
{
    double eps = 1e-6;
    double t = (sqrt(5.0) - 1) / 2;

    double h1 = b - t * (b - a);
    double h2 = a + t * (b - a);

    double f1 = phi(h1, x1, x2, d1, d2);
    double f2 = phi(h2, x1, x2, d1, d2);

    while (b - a > eps)
    {
        if (f1 < f2)
        {
            b = h2;

            h2 = h1;
            f2 = f1;

            h1 = b - t * (b - a);
            f1 = phi(h1, x1, x2, d1, d2);
        }
        else
        {
            a = h1;

            h1 = h2;
            f1 = f2;

            h2 = a + t * (b - a);
            f2 = phi(h2, x1, x2, d1, d2);
        }
    }

    return (a + b) / 2;
}

double lineSearch(double x1, double x2,
                  double d1, double d2)
{
    double a, b;

    svenn(x1, x2, d1, d2, a, b);

    return golden(a, b, x1, x2, d1, d2);
}

int main()
{
    double x1 = 9;
    double x2 = 0;

    double eps = 0.01;

    int k = 0;

    cout << fixed << setprecision(6);

    cout << "Antygradient + tochnyi poshuk\n\n";

    cout << left
         << setw(4)  << "k"
         << setw(11) << "x1"
         << setw(11) << "x2"
         << setw(13) << "F"
         << setw(11) << "g1"
         << setw(11) << "g2"
         << setw(12) << "|g|"
         << setw(12) << "d1"
         << setw(12) << "d2"
         << setw(11) << "h"
         << "Fnew"
         << endl;

    cout << string(115, '-') << endl;

    while (k < 100)
    {
        double g1, g2;
        grad(x1, x2, g1, g2);

        double gn = sqrt(g1*g1 + g2*g2);

        if (gn <= eps)
            break;

        double fx = F(x1, x2);

        // завжди чистий антиградієнт
        double d1 = -g1;
        double d2 = -g2;

        // точний пошук
        double h = lineSearch(x1, x2, d1, d2);

        double xn1 = x1 + h*d1;
        double xn2 = x2 + h*d2;

        double fn = F(xn1, xn2);

        cout << left
             << setw(4)  << k
             << setw(11) << x1
             << setw(11) << x2
             << setw(13) << fx
             << setw(11) << g1
             << setw(11) << g2
             << setw(12) << gn
             << setw(12) << d1
             << setw(12) << d2
             << setw(11) << h
             << fn
             << endl;

        x1 = xn1;
        x2 = xn2;

        k++;
    }

    cout << "\nResult:\n";
    cout << "x* = (" << x1 << "; " << x2 << ")" << endl;
    cout << "F(x*) = " << F(x1, x2) << endl;
    cout << "Iteraciy: " << k << endl;
    cout << "NF = " << nf << endl;
    cout << "NgradF = " << ng << endl;
    cout << "Razom = " << nf + ng << endl;

    return 0;
}