#include <iostream>
#include <iomanip>
#include <cmath>
using namespace std;

int nf = 0;

double F(double x1, double x2)
{
    nf++;

    return (1 - x1) * (1 - x1)
         + 100 * (x2 - x1 * x1) * (x2 - x1 * x1);
}

// дослідний пошук
void explore(double y1, double y2, double Fy, double d,
             double &z1, double &z2, double &Fz, int &pr)
{
    z1 = y1;
    z2 = y2;
    Fz = Fy;

    pr = 0;

    double Ft;

    // x1 +
    Ft = F(z1 + d, z2);
    pr++;

    if (Ft < Fz)
    {
        z1 = z1 + d;
        Fz = Ft;
    }
    else
    {
        // x1 -
        Ft = F(z1 - d, z2);
        pr++;

        if (Ft < Fz)
        {
            z1 = z1 - d;
            Fz = Ft;
        }
    }

    // x2 +
    Ft = F(z1, z2 + d);
    pr++;

    if (Ft < Fz)
    {
        z2 = z2 + d;
        Fz = Ft;
    }
    else
    {
        // x2 -
        Ft = F(z1, z2 - d);
        pr++;

        if (Ft < Fz)
        {
            z2 = z2 - d;
            Fz = Ft;
        }
    }
}

int main()
{
    double x1 = -1.2;
    double x2 = 1.0;

    double d = 1.0;
    double eps = 0.0001;

    double Fb = F(x1, x2);

    double y1 = x1;
    double y2 = x2;
    double Fy = Fb;

    bool pattern = false;

    int k = 0;

    cout << fixed << setprecision(6);

    cout << "Metod Hooka-Jeevesa - Rosenbrock\n\n";

    cout << left
         << setw(4)  << "k"
         << setw(22) << "Diya"
         << setw(12) << "zvidky x1"
         << setw(12) << "zvidky x2"
         << setw(13) << "F"
         << setw(10) << "d"
         << setw(12) << "kudy x1"
         << setw(12) << "kudy x2"
         << setw(13) << "Fnew"
         << setw(8)  << "prob"
         << "NF"
         << endl;

    cout << string(130, '-') << endl;

    while (true)
    {
        double z1, z2, Fz;
        int pr;

        explore(y1, y2, Fy, d, z1, z2, Fz, pr);

        cout << left
             << setw(4)  << k
             << setw(22)
             << (pattern ? "doslidnyi zi zrazka" : "doslidnyi z bazy")
             << setw(12) << y1
             << setw(12) << y2
             << setw(13) << Fy
             << setw(10) << d
             << setw(12) << z1
             << setw(12) << z2
             << setw(13) << Fz
             << setw(8)  << pr
             << nf
             << endl;

        if (Fz < Fb)
        {
            double old1 = x1;
            double old2 = x2;

            x1 = z1;
            x2 = z2;
            Fb = Fz;

            y1 = z1 + (z1 - old1);
            y2 = z2 + (z2 - old2);

            Fy = F(y1, y2);

            k++;

            cout << left
                 << setw(4)  << k
                 << setw(22) << "nova baza -> zrazok"
                 << setw(12) << x1
                 << setw(12) << x2
                 << setw(13) << Fb
                 << setw(10) << d
                 << setw(12) << y1
                 << setw(12) << y2
                 << setw(13) << Fy
                 << setw(8)  << 1
                 << nf
                 << endl;

            pattern = true;
        }
        else if (pattern)
        {
            y1 = x1;
            y2 = x2;
            Fy = Fb;

            pattern = false;

            cout << left
                 << setw(4)  << k
                 << setw(22) << "povernennya v bazu"
                 << setw(12) << x1
                 << setw(12) << x2
                 << setw(13) << Fb
                 << setw(10) << d
                 << setw(12) << "-"
                 << setw(12) << "-"
                 << setw(13) << "-"
                 << setw(8)  << 0
                 << nf
                 << endl;
        }
        else
        {
            d = d / 2;

            cout << left
                 << setw(4)  << k
                 << setw(22) << "droblennya d"
                 << setw(12) << x1
                 << setw(12) << x2
                 << setw(13) << Fb
                 << setw(10) << d
                 << setw(12) << "-"
                 << setw(12) << "-"
                 << setw(13) << "-"
                 << setw(8)  << 0
                 << nf
                 << endl;

            if (d < eps)
                break;

            y1 = x1;
            y2 = x2;
            Fy = Fb;
        }
    }

    cout << "\nResult:\n";
    cout << "x* = (" << x1 << "; " << x2 << ")" << endl;
    cout << "F(x*) = " << Fb << endl;
    cout << "Bazovyh tochok: " << k << endl;
    cout << "NF = " << nf << endl;

    return 0;
}