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
    double eps = 0.001;

    double x1, x2, f1, f2, L;
    int k = 0;
    int n = 0;

    cout << fixed << setprecision(6);

    cout << "Metod dihottomii\n\n";

    cout << left
         << setw(5) << "k"
         << setw(12) << "x1"
         << setw(12) << "x2"
         << setw(15) << "f(x1)"
         << setw(15) << "f(x2)"
         << setw(25) << "[a; b]"
         << setw(12) << "L"
         << endl;

    cout << string(96, '-') << endl;

    L = b - a;

    cout << left
         << setw(5) << 0
         << setw(12) << "-"
         << setw(12) << "-"
         << setw(15) << "-"
         << setw(15) << "-"
         << "[" << a << "; " << b << "]"
         << setw(10) << " "
         << L << endl;

    while (L > sig)
    {
        k++;

        x1 = (a + b) / 2 - eps / 2;
        x2 = x1 + eps;

        f1 = f(x1);
        f2 = f(x2);

        n += 2;

        if (f1 < f2)
        {
            b = x2;
        }
        else
        {
            a = x1;
        }

        L = b - a;

        cout << left
             << setw(5) << k
             << setw(12) << x1
             << setw(12) << x2
             << setw(15) << f1
             << setw(15) << f2
             << "[" << a << "; " << b << "]"
             << setw(10) << " "
             << L << endl;
    }

    double x = (a + b) / 2;

    cout << "\nResult:\n";
    cout << "x* = " << x << endl;
    cout << "f(x*) = " << f(x) << endl;
    cout << "Obchyslen f(x): " << n << endl;

    return 0;
}