#include <iostream>
#include <iomanip>
#include <cmath>
using namespace std;

double f(double x)
{
    return pow(1 - exp(x) * sin(x), 2);
}

int main()
{
    double x0 = 0.75;
    double d = 0.1;

    double f0 = f(x0);
    double x1 = x0 + d;
    double f1 = f(x1);

    if (f1 >= f0)
    {
        d = -d;
        x1 = x0 + d;
        f1 = f(x1);
    }

    cout << fixed << setprecision(6);

    cout << "Metod Svenna\n\n";

    cout << left
         << setw(5) << "k"
         << setw(12) << "d"
         << setw(12) << "xk"
         << setw(15) << "f(xk)"
         << setw(20) << "f(xk)<f(xk-1)"
         << setw(20) << "[a0; b0]"
         << endl;

    cout << string(84, '-') << endl;

    cout << left
         << setw(5) << 0
         << setw(12) << d
         << setw(12) << x0
         << setw(15) << f0
         << setw(20) << "-"
         << setw(20) << "-"
         << endl;

    cout << left
         << setw(5) << 1
         << setw(12) << d
         << setw(12) << x1
         << setw(15) << f1
         << setw(20) << (f1 < f0 ? "+" : "-")
         << setw(20) << "-"
         << endl;

    double xp = x0;
    double x = x1;
    double fx = f1;

    int k = 1;

    while (true)
    {
        k++;

        double xn = x + pow(2, k - 1) * d;
        double fn = f(xn);

        if (fn < fx)
        {
            cout << left
                 << setw(5) << k
                 << setw(12) << d
                 << setw(12) << xn
                 << setw(15) << fn
                 << setw(20) << "+"
                 << setw(20) << "-"
                 << endl;
        }
        else
        {
            double a, b;

            if (d > 0)
            {
                a = xp;
                b = xn;
            }
            else
            {
                a = xn;
                b = xp;
            }

            cout << left
                 << setw(5) << k
                 << setw(12) << d
                 << setw(12) << xn
                 << setw(15) << fn
                 << setw(20) << "-"
                 << "[" << a << "; " << b << "]"
                 << endl;

            cout << "\nL0 = " << b - a << endl;

            break;
        }

        xp = x;
        x = xn;
        fx = fn;
    }

    return 0;
}