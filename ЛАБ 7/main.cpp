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

int main()
{
    double xm = 0.45;
    double xk = 0.75;
    double eps = 0.001;

    double pm = df(xm);
    double pk = df(xk);

    int n = 2;

    cout << fixed << setprecision(6);

    cout << "Metod sichnyh\n\n";

    if (pm * pk >= 0)
    {
        cout << "f' ne zminyue znak na intervali\n";
        return 0;
    }

    cout << left
         << setw(4) << "k"
         << setw(12) << "x(k-1)"
         << setw(15) << "f'(x(k-1))"
         << setw(12) << "xk"
         << setw(15) << "f'(xk)"
         << setw(12) << "x(k+1)"
         << setw(15) << "f'(x(k+1))"
         << setw(13) << "f(x(k+1))"
         << "Diya"
         << endl;

    cout << string(115, '-') << endl;

    double xn = 0;
    double pn = 0;

    for (int k = 1; k <= 100; k++)
    {
        if (abs(pk - pm) < 1e-12)
        {
            cout << "\nDilennya na nul.\n";
            return 0;
        }

        xn = (pk * xm - pm * xk) / (pk - pm);
        pn = df(xn);
        n++;

        cout << left
             << setw(4) << k
             << setw(12) << xm
             << setw(15) << pm
             << setw(12) << xk
             << setw(15) << pk
             << setw(12) << xn
             << setw(15) << pn
             << setw(13) << f(xn);

        if (abs(pn) <= eps)
        {
            cout << "STOP" << endl;
            break;
        }

        if (pn > 0)
        {
            xk = xn;
            pk = pn;
            cout << "xk <- xn" << endl;
        }
        else
        {
            xm = xn;
            pm = pn;
            cout << "xm <- xn" << endl;
        }
    }

    cout << "\nResult:\n";
    cout << "x* = " << xn << endl;

    cout << setprecision(9);
    cout << "f(x*) = " << f(xn) << endl;

    cout << "Obchyslen f'(x): " << n << endl;

    return 0;
}