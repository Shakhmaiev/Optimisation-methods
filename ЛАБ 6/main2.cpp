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
    double alpha = 0.5;
    double tmin = 1e-8;

    int k = 0;

    cout << fixed << setprecision(6);

    cout << "Newton-Rafson z reguluvannyam kroku\n\n";

    cout << left
         << setw(4) << "k"
         << setw(11) << "xk"
         << setw(13) << "f(xk)"
         << setw(13) << "f'(xk)"
         << setw(13) << "f''(xk)"
         << setw(8) << "t"
         << setw(11) << "xb"
         << setw(13) << "f(xb)"
         << setw(13) << "f'(xb)"
         << "Umova"
         << endl;

    cout << string(110, '-') << endl;

    while (k < 50)
    {
        double fx = f(x);
        double d1 = df(x);

        if (abs(d1) <= eps)
            break;

        double d2 = d2f(x);

        if (d2 <= 0)
        {
            cout << "\nf'' <= 0. Zminit x0.\n";
            return 0;
        }

        double t = 1.0;

        while (true)
        {
            double xb = x - t * d1 / d2;
            double d1b = df(xb);

            if (abs(d1b) <= eps)
            {
                cout << left
                     << setw(4) << k + 1
                     << setw(11) << x
                     << setw(13) << fx
                     << setw(13) << d1
                     << setw(13) << d2
                     << setw(8) << t
                     << setw(11) << xb
                     << setw(13) << f(xb)
                     << setw(13) << d1b
                     << "STOP"
                     << endl;

                x = xb;
                k++;
                break;
            }

            double fb = f(xb);

            double leftSide = fb - fx;
            double rightSide =
                -alpha * t * d1 * d1 / d2;

            bool ok = leftSide <= rightSide;

            cout << left
                 << setw(4) << k + 1
                 << setw(11) << x
                 << setw(13) << fx
                 << setw(13) << d1
                 << setw(13) << d2
                 << setw(8) << t
                 << setw(11) << xb
                 << setw(13) << fb
                 << setw(13) << d1b
                 << (ok ? "+" : "-")
                 << endl;

            if (ok)
            {
                x = xb;
                k++;
                break;
            }

            t = t / 2;

            if (t < tmin)
            {
                cout << "\nKrok zanadto malyy.\n";
                return 0;
            }
        }

        if (abs(df(x)) <= eps)
            break;
    }

    cout << "\nResult:\n";
    cout << "x* = " << x << endl;

    cout << setprecision(9);
    cout << "f(x*) = " << f(x) << endl;

    cout << "Iteraciy: " << k << endl;

    return 0;
}