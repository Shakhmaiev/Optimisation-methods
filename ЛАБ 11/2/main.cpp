#include <iostream>
#include <iomanip>
#include <cmath>
using namespace std;

int nPhi = 0;
int nGrad = 0;

double r = 1.0;
double mu = 10.0;
double epsInner = 0.01;

int mode = 1;

/*
mode:
1 - одне обмеження g1 <= 0
2 - рівність h = 0
3 - два обмеження g1 <= 0, g2 <= 0
4 - коло
5 - логарифмічний бар'єр
*/


double F(double x1, double x2)
{
    return 4 * x1 * x1
         + 3 * x2 * x2
         + 2 * x1 * x2
         - 18 * x1
         + 12 * x2;
}

void gradF(double x1, double x2,
           double &g1, double &g2)
{
    g1 = 8 * x1 + 2 * x2 - 18;
    g2 = 2 * x1 + 6 * x2 + 12;
}


double G1(double x1, double x2)
{
    return x1 - x2 - 3;
}

double G2(double x1, double x2)
{
    return x2 + 2;
}

double Gcircle(double x1, double x2)
{
    return 36 - x1 * x1 - x2 * x2;
}


double penaltyValue(double x1, double x2)
{
    if (mode == 1)
    {
        double s = G1(x1, x2);

        if (s < 0)
            s = 0;

        return s * s;
    }

    if (mode == 2)
    {
        double h = G1(x1, x2);

        return h * h;
    }

    if (mode == 3)
    {
        double s1 = G1(x1, x2);
        double s2 = G2(x1, x2);

        if (s1 < 0)
            s1 = 0;

        if (s2 < 0)
            s2 = 0;

        return s1 * s1 + s2 * s2;
    }

    if (mode == 4)
    {
        double s = Gcircle(x1, x2);

        if (s < 0)
            s = 0;

        return s * s;
    }

    return 0;
}


double Phi(double x1, double x2)
{
    nPhi++;

    if (mode == 5)
    {
        double g = G1(x1, x2);

        if (g >= 0)
            return 1e100;

        return F(x1, x2) - mu * log(-g);
    }

    return F(x1, x2)
         + r * penaltyValue(x1, x2);
}


void gradPhi(double x1, double x2,
             double &p1, double &p2)
{
    nGrad++;

    double f1, f2;

    gradF(x1, x2, f1, f2);


    if (mode == 1)
    {
        double s = G1(x1, x2);

        if (s < 0)
            s = 0;


        p1 = f1 + 2 * r * s;
        p2 = f2 - 2 * r * s;

        return;
    }


    if (mode == 2)
    {
        double h = G1(x1, x2);

        p1 = f1 + 2 * r * h;
        p2 = f2 - 2 * r * h;

        return;
    }


    if (mode == 3)
    {
        double s1 = G1(x1, x2);
        double s2 = G2(x1, x2);

        if (s1 < 0)
            s1 = 0;

        if (s2 < 0)
            s2 = 0;


        p1 = f1 + 2 * r * s1;

        p2 = f2
           - 2 * r * s1
           + 2 * r * s2;

        return;
    }


    if (mode == 4)
    {
        double s = Gcircle(x1, x2);

        if (s < 0)
            s = 0;


        p1 = f1 - 4 * r * s * x1;
        p2 = f2 - 4 * r * s * x2;

        return;
    }


    if (mode == 5)
    {
        double g = G1(x1, x2);


        p1 = f1 - mu / g;
        p2 = f2 + mu / g;

        return;
    }
}


double lineF(double x1, double x2,
             double d1, double d2,
             double h)
{
    return Phi(x1 + h * d1,
               x2 + h * d2);
}


void svenn(double x1, double x2,
           double d1, double d2,
           double &a, double &b)
{
    double h0 = 0;
    double delta = 0.01;

    double f0 =
        lineF(x1, x2, d1, d2, h0);

    double f1 =
        lineF(x1, x2, d1, d2, h0 + delta);

    if (f1 > f0)
    {
        delta = -delta;

        f1 =
            lineF(x1, x2,
                  d1, d2,
                  h0 + delta);
    }

    if (f1 >= f0)
    {
        a = h0 + delta;
        b = h0 - delta;

        if (a > b)
        {
            double t = a;
            a = b;
            b = t;
        }

        return;
    }

    double hPrev = h0;
    double h = h0 + delta;

    int k = 1;

    while (true)
    {
        double hNext =
            h + pow(2.0, k) * delta;

        double fh =
            lineF(x1, x2,
                  d1, d2, h);

        double fNext =
            lineF(x1, x2,
                  d1, d2, hNext);

        if (fNext >= fh)
        {
            a = hPrev;
            b = hNext;

            if (a > b)
            {
                double t = a;
                a = b;
                b = t;
            }

            return;
        }

        hPrev = h;
        h = hNext;

        k++;

        if (k > 50)
        {
            a = hPrev;
            b = h;

            return;
        }
    }
}


double golden(double x1, double x2,
              double d1, double d2,
              double a, double b)
{
    double eps = 0.000001;

    double q =
        (sqrt(5.0) - 1) / 2;

    double xL =
        b - q * (b - a);

    double xR =
        a + q * (b - a);

    double fL =
        lineF(x1, x2,
              d1, d2, xL);

    double fR =
        lineF(x1, x2,
              d1, d2, xR);

    while (fabs(b - a) > eps)
    {
        if (fL < fR)
        {
            b = xR;

            xR = xL;
            fR = fL;

            xL =
                b - q * (b - a);

            fL =
                lineF(x1, x2,
                      d1, d2, xL);
        }
        else
        {
            a = xL;

            xL = xR;
            fL = fR;

            xR =
                a + q * (b - a);

            fR =
                lineF(x1, x2,
                      d1, d2, xR);
        }
    }

    return (a + b) / 2;
}


int fletcherReeves(double &x1,
                   double &x2)
{
    double g1, g2;

    gradPhi(x1, x2, g1, g2);

    double d1 = -g1;
    double d2 = -g2;

    int iter = 0;

    while (true)
    {
        double norm =
            sqrt(g1 * g1 + g2 * g2);

        if (norm <= epsInner)
            break;

        double a, b;

        svenn(x1, x2,
              d1, d2,
              a, b);

        double h =
            golden(x1, x2,
                   d1, d2,
                   a, b);

        x1 += h * d1;
        x2 += h * d2;

        double ng1, ng2;

        gradPhi(x1, x2,
                ng1, ng2);

        double beta =
            (ng1 * ng1 + ng2 * ng2) /
            (g1 * g1 + g2 * g2);

        iter++;

        if (iter % 2 == 0)
            beta = 0;

        d1 = -ng1 + beta * d1;
        d2 = -ng2 + beta * d2;

        g1 = ng1;
        g2 = ng2;

        if (iter >= 5000)
            break;
    }

    return iter;
}


struct Result
{
    double x1;
    double x2;

    int outer;
    int totalPhi;
    int totalGrad;

    double violation;
    double f;
    double lambda1;
    double lambda2;
};


Result SUMT(double start1,
            double start2,
            double c,
            double innerEps,
            bool printTable)
{
    double x1 = start1;
    double x2 = start2;

    r = 1.0;
    epsInner = innerEps;

    int outer = 0;
    int totalPhi = 0;
    int totalGrad = 0;

    if (printTable)
    {
        cout << left
             << setw(5)  << "k"
             << setw(12) << "r"
             << setw(8)  << "iter"
             << setw(10) << "NPhi"
             << setw(10) << "Ngrad"
             << setw(13) << "x1"
             << setw(13) << "x2"
             << setw(14) << "F"
             << setw(14) << "porush"
             << setw(14) << "rP"
             << endl;
    }

    while (true)
    {
        nPhi = 0;
        nGrad = 0;

        int it =
            fletcherReeves(x1, x2);

        double P =
            penaltyValue(x1, x2);

        double rP = r * P;

        double violation = 0;

        if (mode == 1)
            violation =
                max(0.0, G1(x1, x2));

        if (mode == 2)
            violation =
                fabs(G1(x1, x2));

        if (mode == 3)
            violation =
                max(max(0.0, G1(x1, x2)),
                    max(0.0, G2(x1, x2)));

        if (mode == 4)
            violation =
                max(0.0, Gcircle(x1, x2));

        totalPhi += nPhi;
        totalGrad += nGrad;

        outer++;

        if (printTable)
        {
            cout << left
                 << setw(5)  << outer
                 << setw(12) << r
                 << setw(8)  << it
                 << setw(10) << nPhi
                 << setw(10) << nGrad
                 << setw(13) << x1
                 << setw(13) << x2
                 << setw(14) << F(x1, x2)
                 << setw(14) << violation
                 << setw(14) << rP
                 << endl;
        }

        if (rP <= 0.01)
            break;

        r *= c;

        if (outer > 100)
            break;
    }

    Result res;

    res.x1 = x1;
    res.x2 = x2;

    res.outer = outer;
    res.totalPhi = totalPhi;
    res.totalGrad = totalGrad;

    res.f = F(x1, x2);

    res.lambda1 = 0;
    res.lambda2 = 0;

    if (mode == 1)
    {
        res.violation =
            max(0.0, G1(x1, x2));

        res.lambda1 =
            2 * r *
            max(0.0, G1(x1, x2));
    }

    if (mode == 2)
    {
        res.violation =
            fabs(G1(x1, x2));

        res.lambda1 =
            2 * r * G1(x1, x2);
    }

    if (mode == 3)
    {
        double s1 =
            max(0.0, G1(x1, x2));

        double s2 =
            max(0.0, G2(x1, x2));

        res.violation =
            max(s1, s2);

        res.lambda1 = 2 * r * s1;
        res.lambda2 = 2 * r * s2;
    }

    if (mode == 4)
    {
        res.violation =
            max(0.0,
                Gcircle(x1, x2));

        res.lambda1 =
            2 * r *
            max(0.0,
                Gcircle(x1, x2));
    }

    return res;
}


void barrierMethod()
{
    mode = 5;

    double x1 = 0;
    double x2 = 0;

    mu = 10.0;
    epsInner = 0.01;

    int totalPhi = 0;
    int totalGrad = 0;

    cout << left
         << setw(12) << "mu"
         << setw(13) << "x1"
         << setw(13) << "x2"
         << setw(14) << "F"
         << setw(14) << "g"
         << setw(14) << "lambda"
         << setw(10) << "NPhi"
         << setw(10) << "Ngrad"
         << endl;

    while (true)
    {
        nPhi = 0;
        nGrad = 0;

        fletcherReeves(x1, x2);

        totalPhi += nPhi;
        totalGrad += nGrad;

        double g =
            G1(x1, x2);

        double lambda =
            mu / (-g);

        cout << left
             << setw(12) << mu
             << setw(13) << x1
             << setw(13) << x2
             << setw(14) << F(x1, x2)
             << setw(14) << g
             << setw(14) << lambda
             << setw(10) << nPhi
             << setw(10) << nGrad
             << endl;

        if (mu <= 0.01)
            break;

        mu /= 10.0;
    }

    cout << "NPhi sum = "
         << totalPhi << endl;

    cout << "Ngrad sum = "
         << totalGrad << endl;
}


int sphereMethod(double &x1,
                 double &x2)
{
    double h = 1.0;

    int iter = 0;

    double Fx =
        Phi(x1, x2);

    while (true)
    {
        double g1, g2;

        gradPhi(x1, x2,
                g1, g2);

        double norm =
            sqrt(g1 * g1 +
                 g2 * g2);

        if (norm <= 0.01)
            break;

        double d1 = -g1 / norm;
        double d2 = -g2 / norm;

        double nx1 =
            x1 + h * d1;

        double nx2 =
            x2 + h * d2;

        double Fn =
            Phi(nx1, nx2);

        while (Fn >= Fx)
        {
            h /= 2.0;

            if (h < 1e-10)
                return iter;

            nx1 =
                x1 + h * d1;

            nx2 =
                x2 + h * d2;

            Fn =
                Phi(nx1, nx2);
        }

        x1 = nx1;
        x2 = nx2;
        Fx = Fn;

        iter++;

        if (iter >= 100000)
            break;
    }

    return iter;
}


int main()
{
    cout << fixed
         << setprecision(6);


    cout << "\n------------------------------\n";
    cout << "1. EXPERIMENT C\n";
    cout << "------------------------------\n";

    mode = 1;

    double cs[4] =
    {
        2, 4, 10, 100
    };

    cout << left
         << setw(10) << "c"
         << setw(12) << "outer"
         << setw(12) << "NPhi"
         << setw(12) << "Ngrad"
         << setw(13) << "x1"
         << setw(13) << "x2"
         << setw(14) << "F"
         << setw(14) << "porush"
         << endl;

    for (int i = 0; i < 4; i++)
    {
        Result a =
            SUMT(9, 0,
                 cs[i],
                 0.01,
                 false);

        cout << left
             << setw(10) << cs[i]
             << setw(12) << a.outer
             << setw(12) << a.totalPhi
             << setw(12) << a.totalGrad
             << setw(13) << a.x1
             << setw(13) << a.x2
             << setw(14) << a.f
             << setw(14) << a.violation
             << endl;
    }


    cout << "\n------------------------------\n";
    cout << "2. EXPERIMENT EPS\n";
    cout << "------------------------------\n";

    double epsArray[3] =
    {
        0.1,
        0.01,
        0.001
    };

    cout << left
         << setw(12) << "eps"
         << setw(12) << "outer"
         << setw(12) << "NPhi"
         << setw(12) << "Ngrad"
         << setw(13) << "x1"
         << setw(13) << "x2"
         << setw(14) << "F"
         << setw(14) << "porush"
         << endl;

    for (int i = 0; i < 3; i++)
    {
        mode = 1;

        Result a =
            SUMT(9, 0,
                 10,
                 epsArray[i],
                 false);

        cout << left
             << setw(12) << epsArray[i]
             << setw(12) << a.outer
             << setw(12) << a.totalPhi
             << setw(12) << a.totalGrad
             << setw(13) << a.x1
             << setw(13) << a.x2
             << setw(14) << a.f
             << setw(14) << a.violation
             << endl;
    }


    cout << "\n------------------------------\n";
    cout << "3. RIVNIST x1 - x2 = 3\n";
    cout << "------------------------------\n";

    mode = 2;

    Result eq =
        SUMT(9, 0,
             10,
             0.01,
             true);

    cout << "\nResult:\n";

    cout << "x = ("
         << eq.x1 << "; "
         << eq.x2 << ")\n";

    cout << "F = "
         << eq.f << endl;

    cout << "|h| = "
         << eq.violation << endl;

    cout << "lambda = "
         << eq.lambda1 << endl;

    cout << "NPhi = "
         << eq.totalPhi << endl;


    cout << "\n------------------------------\n";
    cout << "4. LOG-BARRIER\n";
    cout << "------------------------------\n";

    barrierMethod();


    cout << "\n------------------------------\n";
    cout << "5. DVA OBMEZHENNYA\n";
    cout << "g1 = x1-x2-3 <= 0\n";
    cout << "g2 = x2+2 <= 0\n";
    cout << "------------------------------\n";

    mode = 3;

    Result two =
        SUMT(9, 0,
             10,
             0.01,
             true);

    cout << "\nResult:\n";

    cout << "x = ("
         << two.x1 << "; "
         << two.x2 << ")\n";

    cout << "F = "
         << two.f << endl;

    cout << "g1 = "
         << G1(two.x1,
               two.x2)
         << endl;

    cout << "g2 = "
         << G2(two.x1,
               two.x2)
         << endl;

    cout << "lambda1 = "
         << two.lambda1
         << endl;

    cout << "lambda2 = "
         << two.lambda2
         << endl;

    cout << "NPhi = "
         << two.totalPhi
         << endl;


    cout << "\n------------------------------\n";
    cout << "6. KOLO\n";
    cout << "x1^2 + x2^2 >= 36\n";
    cout << "------------------------------\n";

    mode = 4;

    cout << "\nStart 1: (9; 0)\n";

    Result circle1 =
        SUMT(9, 0,
             10,
             0.01,
             false);

    cout << "x = ("
         << circle1.x1 << "; "
         << circle1.x2 << ")\n";

    cout << "F = "
         << circle1.f << endl;

    cout << "norm = "
         << sqrt(circle1.x1 *
                 circle1.x1 +
                 circle1.x2 *
                 circle1.x2)
         << endl;

    cout << "porush = "
         << circle1.violation
         << endl;

    cout << "NPhi = "
         << circle1.totalPhi
         << endl;


    cout << "\nStart 2: (-9; 0)\n";

    Result circle2 =
        SUMT(-9, 0,
             10,
             0.01,
             false);

    cout << "x = ("
         << circle2.x1 << "; "
         << circle2.x2 << ")\n";

    cout << "F = "
         << circle2.f << endl;

    cout << "norm = "
         << sqrt(circle2.x1 *
                 circle2.x1 +
                 circle2.x2 *
                 circle2.x2)
         << endl;

    cout << "porush = "
         << circle2.violation
         << endl;

    cout << "NPhi = "
         << circle2.totalPhi
         << endl;


    cout << "\n------------------------------\n";
    cout << "7. FR VS SPHERE, r = 100\n";
    cout << "------------------------------\n";

    mode = 1;
    r = 100;
    epsInner = 0.01;


    double fx1 = 9;
    double fx2 = 0;

    nPhi = 0;
    nGrad = 0;

    int frIter =
        fletcherReeves(fx1,
                       fx2);

    int frPhi = nPhi;
    int frGrad = nGrad;

    cout << "\nFletcher-Reeves:\n";

    cout << "x = ("
         << fx1 << "; "
         << fx2 << ")\n";

    cout << "F = "
         << F(fx1, fx2)
         << endl;

    cout << "iter = "
         << frIter << endl;

    cout << "NPhi = "
         << frPhi << endl;

    cout << "Ngrad = "
         << frGrad << endl;



    double sx1 = 9;
    double sx2 = 0;

    nPhi = 0;
    nGrad = 0;

    int sphereIter =
        sphereMethod(sx1,
                     sx2);

    cout << "\nSpherical metric:\n";

    cout << "x = ("
         << sx1 << "; "
         << sx2 << ")\n";

    cout << "F = "
         << F(sx1, sx2)
         << endl;

    cout << "iter = "
         << sphereIter
         << endl;

    cout << "NPhi = "
         << nPhi << endl;

    cout << "Ngrad = "
         << nGrad << endl;


    cout << "\n------------------------------\n";
    cout << "8. SUMMARY\n";
    cout << "------------------------------\n";

    cout << "Equality:\n";
    cout << "x=("
         << eq.x1 << ";"
         << eq.x2
         << "), F="
         << eq.f
         << ", NPhi="
         << eq.totalPhi
         << endl;

    cout << "\nTwo constraints:\n";
    cout << "x=("
         << two.x1 << ";"
         << two.x2
         << "), F="
         << two.f
         << ", lambda=("
         << two.lambda1 << ";"
         << two.lambda2
         << "), NPhi="
         << two.totalPhi
         << endl;

    cout << "\nCircle start 1:\n";
    cout << "x=("
         << circle1.x1 << ";"
         << circle1.x2
         << "), F="
         << circle1.f
         << ", NPhi="
         << circle1.totalPhi
         << endl;

    cout << "\nCircle start 2:\n";
    cout << "x=("
         << circle2.x1 << ";"
         << circle2.x2
         << "), F="
         << circle2.f
         << ", NPhi="
         << circle2.totalPhi
         << endl;

    cout << "\nFR r=100:\n";
    cout << "iter="
         << frIter
         << ", NPhi="
         << frPhi
         << endl;

    cout << "\nSphere r=100:\n";
    cout << "iter="
         << sphereIter
         << ", NPhi="
         << nPhi
         << endl;

    return 0;
}