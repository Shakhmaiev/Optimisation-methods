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
    double eps = 0.01;

    double t = (sqrt(5.0) - 1) / 2;

    double x1 = a + (1 - t) * (b - a);
    double x2 = a + t * (b - a);

    double f1 = f(x1);
    double f2 = f(x2);

    double L;
    int k = 1;
    int n = 2;

    cout << fixed << setprecision(6);

    cout << "Metod zolotogo pererizu\n\n";

    cout << left
         << setw(5) << "k"
         << setw(12) << "a"
         << setw(12) << "b"
         << setw(12) << "x1"
         << setw(12) << "x2"
         << setw(15) << "f(x1)"
         << setw(15) << "f(x2)"
         << setw(12) << "L"
         << setw(15) << "Rishennya"
         << endl;

    cout << string(110, '-') << endl;

    while (true)
    {
        L = b - a;

        cout << left
             << setw(5) << k
             << setw(12) << a
             << setw(12) << b
             << setw(12) << x1
             << setw(12) << x2
             << setw(15) << f1
             << setw(15) << f2
             << setw(12) << L;

        if (L <= eps)
        {
            cout << "STOP" << endl;
            break;
        }

        if (f1 <= f2)
        {
            cout << "f1 <= f2" << endl;

            b = x2;

            x2 = x1;
            f2 = f1;

            x1 = a + (1 - t) * (b - a);
            f1 = f(x1);

            n++;
        }
        else
        {
            cout << "f1 > f2" << endl;

            a = x1;

            x1 = x2;
            f1 = f2;

            x2 = a + t * (b - a);
            f2 = f(x2);

            n++;
        }

        k++;
    }

    double x = (a + b) / 2;

    cout << "\nResult:\n";
    cout << "x* = " << x << endl;

    cout << setprecision(9);
    cout << "f(x*) = " << f(x) << endl;

    cout << "Obchyslen f(x): " << n << endl;

    return 0;
}