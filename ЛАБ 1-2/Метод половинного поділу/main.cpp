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
    double a = 0.45;
    double b = 0.75;
    double sig = 0.01;

    double L = b - a;

    double xm = (a + b) / 2;
    double fm = f(xm);

    double x1, x2, f1, f2;

    int k = 0;
    int n = 1;

    cout << fixed << setprecision(6);

    cout << "Metod polovynnogo podilu\n\n";

    cout << left
         << setw(5) << "k"
         << setw(12) << "x1"
         << setw(12) << "xm"
         << setw(12) << "x2"
         << setw(15) << "f(x1)"
         << setw(15) << "f(xm)"
         << setw(15) << "f(x2)"
         << setw(12) << "L"
         << endl;

    cout << string(98, '-') << endl;

    cout << left
         << setw(5) << 0
         << setw(12) << "-"
         << setw(12) << "-"
         << setw(12) << "-"
         << setw(15) << "-"
         << setw(15) << "-"
         << setw(15) << "-"
         << setw(12) << L
         << endl;

    while (L > sig)
    {
        k++;

        x1 = a + L / 4;
        x2 = b - L / 4;

        f1 = f(x1);
        f2 = f(x2);

        n += 2;

        cout << left
             << setw(5) << k
             << setw(12) << x1
             << setw(12) << xm
             << setw(12) << x2
             << setw(15) << f1
             << setw(15) << fm
             << setw(15) << f2;

        if (f1 < fm)
        {
            b = xm;
            xm = x1;
            fm = f1;
        }
        else if (f2 < fm)
        {
            a = xm;
            xm = x2;
            fm = f2;
        }
        else
        {
            a = x1;
            b = x2;
        }

        L = b - a;

        cout << setw(12) << L << endl;
    }

    cout << "\nResult:\n";
    cout << "x* = " << xm << endl;

    cout << setprecision(9);
    cout << "f(x*) = " << fm << endl;

    cout << "Obchyslen f(x): " << n << endl;

    return 0;
}