#include <iostream>
#include <iomanip>
#include <cmath>
using namespace std;

int nPhi = 0;
int nGrad = 0;

double r;

// Цільова функція, варіант 15
double F(double x1, double x2)
{
    return 4 * x1 * x1
         + 3 * x2 * x2
         + 2 * x1 * x2
         - 18 * x1
         + 12 * x2;
}

// Обмеження:
// x1 - x2 <= 3
// g(x) = x1 - x2 - 3 <= 0
double G(double x1, double x2)
{
    return x1 - x2 - 3;
}

// Штрафна функція
double Phi(double x1, double x2)
{
    nPhi++;

    double g = G(x1, x2);
    double s = 0;

    if (g > 0)
        s = g;

    return F(x1, x2) + r * s * s;
}

// Градієнт штрафної функції
void gradPhi(double x1, double x2, double &g1, double &g2)
{
    nGrad++;

    double s = G(x1, x2);

    if (s < 0)
        s = 0;

    g1 = 8 * x1 + 2 * x2 - 18 + 2 * r * s;
    g2 = 2 * x1 + 6 * x2 + 12 - 2 * r * s;
}

// Функція вздовж напрямку
double lineF(double x1, double x2,
             double d1, double d2,
             double h)
{
    return Phi(x1 + h * d1,
               x2 + h * d2);
}

// Метод Свенна
void svenn(double x1, double x2,
           double d1, double d2,
           double &a, double &b)
{
    double h0 = 0;
    double delta = 0.01;

    double f0 = lineF(x1, x2, d1, d2, h0);
    double f1 = lineF(x1, x2, d1, d2, h0 + delta);

    if (f1 > f0)
    {
        delta = -delta;
        f1 = lineF(x1, x2, d1, d2, h0 + delta);
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
        double hNext = h + pow(2.0, k) * delta;

        double fh = lineF(x1, x2, d1, d2, h);
        double fNext = lineF(x1, x2, d1, d2, hNext);

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
    }
}

// Золотий переріз
double golden(double x1, double x2,
              double d1, double d2,
              double a, double b)
{
    double eps = 0.000001;

    double k = (sqrt(5.0) - 1) / 2;

    double xL = b - k * (b - a);
    double xR = a + k * (b - a);

    double fL = lineF(x1, x2, d1, d2, xL);
    double fR = lineF(x1, x2, d1, d2, xR);

    while (fabs(b - a) > eps)
    {
        if (fL < fR)
        {
            b = xR;

            xR = xL;
            fR = fL;

            xL = b - k * (b - a);
            fL = lineF(x1, x2, d1, d2, xL);
        }
        else
        {
            a = xL;

            xL = xR;
            fL = fR;

            xR = a + k * (b - a);
            fR = lineF(x1, x2, d1, d2, xR);
        }
    }

    return (a + b) / 2;
}

// Метод Флетчера-Рівса
int fletcherReeves(double &x1, double &x2)
{
    double eps = 0.01;

    double g1, g2;

    gradPhi(x1, x2, g1, g2);

    double d1 = -g1;
    double d2 = -g2;

    int iter = 0;

    while (true)
    {
        double norm = sqrt(g1 * g1 + g2 * g2);

        if (norm <= eps)
            break;

        double a, b;

        svenn(x1, x2, d1, d2, a, b);

        double h = golden(x1, x2,
                          d1, d2,
                          a, b);

        x1 = x1 + h * d1;
        x2 = x2 + h * d2;

        double ng1, ng2;

        gradPhi(x1, x2, ng1, ng2);

        double beta =
            (ng1 * ng1 + ng2 * ng2) /
            (g1 * g1 + g2 * g2);

        iter++;

        // Рестарт через n = 2 кроки
        if (iter % 2 == 0)
            beta = 0;

        d1 = -ng1 + beta * d1;
        d2 = -ng2 + beta * d2;

        g1 = ng1;
        g2 = ng2;

        // Захист від нескінченного циклу
        if (iter > 1000)
            break;
    }

    return iter;
}

int main()
{
    cout << fixed << setprecision(6);

    // -------------------------------------------------
    // 1. Основний метод зовнішніх штрафних функцій
    // -------------------------------------------------

    double x1 = 9;
    double x2 = 0;

    r = 1.0;

    double c = 10.0;
    double eps = 0.01;

    int totalPhi = 0;
    int totalGrad = 0;

    cout << "Metod zovnishnih shtrafnih funkciy\n\n";

    cout << left
         << setw(4)  << "k"
         << setw(12) << "r"
         << setw(13) << "start x1"
         << setw(13) << "start x2"
         << setw(8)  << "iter"
         << setw(10) << "NPhi"
         << setw(10) << "Ngrad"
         << setw(13) << "x1"
         << setw(13) << "x2"
         << setw(14) << "F"
         << setw(14) << "g"
         << setw(14) << "rP"
         << setw(14) << "lambda"
         << "NPhi sum"
         << endl;

    cout << string(160, '-') << endl;

    int k = 1;

    while (true)
    {
        double start1 = x1;
        double start2 = x2;

        nPhi = 0;
        nGrad = 0;

        int iter = fletcherReeves(x1, x2);

        double g = G(x1, x2);

        double s = 0;

        if (g > 0)
            s = g;

        double P = s * s;
        double rP = r * P;
        double lambda = 2 * r * s;

        totalPhi += nPhi;
        totalGrad += nGrad;

        cout << left
             << setw(4)  << k
             << setw(12) << r
             << setw(13) << start1
             << setw(13) << start2
             << setw(8)  << iter
             << setw(10) << nPhi
             << setw(10) << nGrad
             << setw(13) << x1
             << setw(13) << x2
             << setw(14) << F(x1, x2)
             << setw(14) << g
             << setw(14) << rP
             << setw(14) << lambda
             << totalPhi
             << endl;

        if (rP <= eps)
            break;

        r = r * c;
        k++;
    }

    cout << "\nResult:\n";

    cout << "x = (" << x1
         << "; " << x2 << ")" << endl;

    cout << "F = "
         << F(x1, x2) << endl;

    cout << "g = "
         << G(x1, x2) << endl;

    cout << "NPhi = "
         << totalPhi << endl;

    cout << "NgradPhi = "
         << totalGrad << endl;

    cout << "\nTochna perevirka:\n";

    cout << "x* = (1.666667; -1.333333)"
         << endl;

    cout << "F* = -34.000000"
         << endl;

    cout << "lambda* = 7.333333"
         << endl;

    // -------------------------------------------------
    // 2. Дошліфовування проєкцією
    // -------------------------------------------------

    double gv = G(x1, x2);

    double xp1 = x1 - gv / 2.0;
    double xp2 = x2 + gv / 2.0;

    cout << "\nPislya proekcii:\n";

    cout << "x = (" << xp1
         << "; " << xp2 << ")" << endl;

    cout << "F = "
         << F(xp1, xp2) << endl;

    cout << "g = "
         << G(xp1, xp2) << endl;

    // -------------------------------------------------
    // 3. Один запуск одразу при r = 10^6
    // -------------------------------------------------

    cout << "\n----------------------------------\n";
    cout << "Odyn zapusk pri r = 1000000\n\n";

    x1 = 9;
    x2 = 0;

    r = 1000000.0;

    nPhi = 0;
    nGrad = 0;

    int iter = fletcherReeves(x1, x2);

    double g = G(x1, x2);

    double s = 0;

    if (g > 0)
        s = g;

    double P = s * s;
    double rP = r * P;
    double lambda = 2 * r * s;

    cout << "Result:\n";

    cout << "x = (" << x1
         << "; " << x2 << ")" << endl;

    cout << "F = "
         << F(x1, x2) << endl;

    cout << "g = "
         << g << endl;

    cout << "rP = "
         << rP << endl;

    cout << "lambda = "
         << lambda << endl;

    cout << "Iteraciy = "
         << iter << endl;

    cout << "NPhi = "
         << nPhi << endl;

    cout << "NgradPhi = "
         << nGrad << endl;

    return 0;
}