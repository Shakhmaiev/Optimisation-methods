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
    double x1 = 0.45;
    double dx = 0.1;

    double sf = 0.01;
    double sx = 0.001;

    double x2 = x1 + dx;

    double f1 = f(x1);
    double f2 = f(x2);

    double x3;

    if (f1 > f2)
        x3 = x1 + 2 * dx;
    else
        x3 = x1 - dx;

    double f3 = f(x3);

    int nf = 3;

    cout << fixed << setprecision(6);

    cout << "Metod Pauella\n\n";

    cout << left
         << setw(4) << "k"
         << setw(11) << "x1"
         << setw(11) << "x2"
         << setw(11) << "x3"
         << setw(13) << "f(x1)"
         << setw(13) << "f(x2)"
         << setw(13) << "f(x3)"
         << setw(12) << "a1"
         << setw(12) << "a2"
         << setw(11) << "xt"
         << setw(13) << "f(xt)"
         << setw(12) << "df"
         << setw(12) << "dx"
         << "Diya" << endl;

    cout << string(155, '-') << endl;

    double xt = 0, ft = 0;

    for (int k = 1; k <= 50; k++)
    {
        // перевірка, щоб точки не співпали
        if (abs(x2 - x1) < 1e-12 ||
            abs(x3 - x1) < 1e-12 ||
            abs(x3 - x2) < 1e-12)
        {
            cout << "\nTochky spivpaly. Zminit dx.\n";
            return 0;
        }

       
        double xmin = x1;
        double fmin = f1;

        if (f2 < fmin)
        {
            xmin = x2;
            fmin = f2;
        }

        if (f3 < fmin)
        {
            xmin = x3;
            fmin = f3;
        }

        double xmax = x1;
        double fmax = f1;
        int imax = 1;

        if (f2 > fmax)
        {
            xmax = x2;
            fmax = f2;
            imax = 2;
        }

        if (f3 > fmax)
        {
            xmax = x3;
            fmax = f3;
            imax = 3;
        }

        
        double a1 = (f2 - f1) / (x2 - x1);

        double a2 = ((f3 - f1) / (x3 - x1) - a1)
                    / (x3 - x2);

        
        if (abs(a2) < 1e-12)
        {
            cout << "\na2 = 0. Zminit dx.\n";
            return 0;
        }

        
        if (a2 > 0)
            xt = -a1 / (2 * a2) + (x1 + x2) / 2;
        else
            xt = 2 * xmin - xmax;

        ft = f(xt);
        nf++;

        double df = abs(fmin - ft);
        double dxx = abs(xmin - xt);

        cout << left
             << setw(4) << k
             << setw(11) << x1
             << setw(11) << x2
             << setw(11) << x3
             << setw(13) << f1
             << setw(13) << f2
             << setw(13) << f3
             << setw(12) << a1
             << setw(12) << a2
             << setw(11) << xt
             << setw(13) << ft
             << setw(12) << df
             << setw(12) << dxx;

        
        if (df <= sf && dxx <= sx)
        {
            cout << "STOP" << endl;
            break;
        }

        // якщо нова точка співпала зі старою
        if (abs(xt - x1) < 1e-12 ||
            abs(xt - x2) < 1e-12 ||
            abs(xt - x3) < 1e-12)
        {
            cout << "xt spivpala z tochkoyu -> STOP" << endl;
            break;
        }

        
        if (imax == 1)
        {
            x1 = xt;
            f1 = ft;
            cout << "x1 <- xt" << endl;
        }
        else if (imax == 2)
        {
            x2 = xt;
            f2 = ft;
            cout << "x2 <- xt" << endl;
        }
        else
        {
            x3 = xt;
            f3 = ft;
            cout << "x3 <- xt" << endl;
        }
    }

    cout << "\nResult:\n";
    cout << "x* = " << xt << endl;

    cout << setprecision(9);
    cout << "f(x*) = " << ft << endl;

    cout << "Obchyslen f(x): " << nf << endl;

    return 0;
}