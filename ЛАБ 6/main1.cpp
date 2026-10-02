#include <iostream>
#include <iomanip>
#include <cmath>
using namespace std;

double f(double x)
{
    return pow(1 - exp(x) * sin(x), 2);
}

double df(double x)
{
    return -2 * exp(x) * (sin(x) + cos(x))
           * (1 - exp(x) * sin(x));
}

double d2f(double x)
{
    return 2 * exp(2 * x) * pow(sin(x) + cos(x), 2)
           - 4 * exp(x) * cos(x)
           * (1 - exp(x) * sin(x));
}

int main()
{
    double x = 0.45;
    double eps = 0.001;

    int k = 0;
    int n1 = 0;
    int n2 = 0;

    cout << fixed << setprecision(6);

    cout << "Metod Newtona\n\n";

    cout << left
         << setw(5) << "k"
         << setw(12) << "xk"
         << setw(15) << "f(xk)"
         << setw(15) << "f'(xk)"
         << setw(15) << "|f'|"
         << setw(15) << "f''(xk)"
         << endl;

    cout << string(77, '-') << endl;

    while (true)
    {
        double fx = f(x);
        double d1 = df(x);
        n1++;

        cout << left
             << setw(5) << k
             << setw(12) << x
             << setw(15) << fx
             << setw(15) << d1
             << setw(15) << abs(d1);

        if (abs(d1) <= eps)
        {
            cout << setw(15) << "-" << endl;
            break;
        }

        double d2 = d2f(x);
        n2++;

        cout << setw(15) << d2 << endl;

        if (d2 <= 0)
        {
            cout << "\nf'' <= 0. Zminit x0.\n";
            return 0;
        }

        x = x - d1 / d2;
        k++;
    }

    cout << "\nResult:\n";
    cout << "x* = " << x << endl;
    cout << setprecision(9);
    cout << "f(x*) = " << f(x) << endl;

    cout << "Iteraciy: " << k << endl;
    cout << "f' obchyslen: " << n1 << endl;
    cout << "f'' obchyslen: " << n2 << endl;

    return 0;
}