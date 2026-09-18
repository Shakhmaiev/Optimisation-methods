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
    double d = eps / 10;

    double F[50];
    F[1] = 1;
    F[2] = 1;

    double R = (b - a) / eps;

    int N = 2;

    while (F[N] < R)
    {
        N++;
        F[N] = F[N - 1] + F[N - 2];
    }

    cout << fixed << setprecision(6);

    cout << "Metod Fibonachchi\n\n";
    cout << "R = " << R << endl;
    cout << "N = " << N << endl;
    cout << "F[N] = " << F[N] << endl;
    cout << "delta = " << d << "\n\n";

    double x1 = a + F[N - 2] / F[N] * (b - a);
    double x2 = a + F[N - 1] / F[N] * (b - a);

    double f1 = f(x1);
    double f2 = f(x2);

    int n = 2;

    cout << left
         << setw(4) << "k"
         << setw(4) << "m"
         << setw(12) << "F1/Fm"
         << setw(12) << "F2/Fm"
         << setw(11) << "a"
         << setw(11) << "b"
         << setw(11) << "x1"
         << setw(11) << "x2"
         << setw(13) << "f(x1)"
         << setw(13) << "f(x2)"
         << setw(11) << "L"
         << "Rishennya" << endl;

    cout << string(125, '-') << endl;

    for (int k = 0; k <= N - 3; k++)
    {
        int m = N - k;
        double L = b - a;

        cout << left
             << setw(4) << k + 1
             << setw(4) << m;

        if (m > 3)
        {
            cout << setw(12) << F[m - 2] / F[m]
                 << setw(12) << F[m - 1] / F[m];
        }
        else
        {
            cout << setw(12) << 0.5
                 << setw(12) << 0.5;
        }

        cout << setw(11) << a
             << setw(11) << b
             << setw(11) << x1
             << setw(11) << x2
             << setw(13) << f1
             << setw(13) << f2
             << setw(11) << L;

        if (f1 <= f2)
        {
            cout << "f1 <= f2" << endl;

            b = x2;

            if (m > 3)
            {
                x2 = x1;
                f2 = f1;

                if (m - 1 == 3)
                    x1 = x2 - d;
                else
                    x1 = a + F[m - 3] / F[m - 1] * (b - a);

                f1 = f(x1);
                n++;
            }
        }
        else
        {
            cout << "f1 > f2" << endl;

            a = x1;

            if (m > 3)
            {
                x1 = x2;
                f1 = f2;

                if (m - 1 == 3)
                    x2 = x1 + d;
                else
                    x2 = a + F[m - 2] / F[m - 1] * (b - a);

                f2 = f(x2);
                n++;
            }
        }
    }

    double x = (a + b) / 2;

    cout << "\nResult:\n";
    cout << "Interval: [" << a << "; " << b << "]" << endl;
    cout << "L = " << b - a << endl;
    cout << "x* = " << x << endl;

    cout << setprecision(9);
    cout << "f(x*) = " << f(x) << endl;

    cout << "Obchyslen f(x): " << n << endl;

    return 0;
}