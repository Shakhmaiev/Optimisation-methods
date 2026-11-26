#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <limits>
#include <string>

using namespace std;

struct Wall
{
    double a1, a2, b;
    string name;
};

struct Result
{
    double x1, x2;

    vector<int> active;
    vector<double> u;

    int steps;
    int NF;
    int Ngrad;
};

int NF = 0;
int Ngrad = 0;

const double EPS = 0.01;
const double ACTIVE_EPS = 1e-8;
const double INF = 1e100;


// Цільова функція, варіант 15


double F(double x1, double x2)
{
    NF++;

    return 4 * x1 * x1
         + 3 * x2 * x2
         + 2 * x1 * x2
         - 18 * x1
         + 12 * x2;
}

double Fraw(double x1, double x2)
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
    Ngrad++;

    g1 = 8 * x1 + 2 * x2 - 18;
    g2 = 2 * x1 + 6 * x2 + 12;
}

 
// g(x) = a^T x  b


double G(const Wall &w,
         double x1, double x2)
{
    return w.a1 * x1
         + w.a2 * x2
         - w.b;
}

bool feasible(double x1, double x2,
              const vector<Wall> &walls,
              double eps = 1e-9)
{
    for (const Wall &w : walls)
    {
        if (G(w, x1, x2) > eps)
            return false;
    }

    return true;
}


// Проєкція початкової точки на допустиму множину


void projectPoint(double &x1, double &x2,
                  const vector<Wall> &walls)
{
    if (feasible(x1, x2, walls))
        return;

    double ox = x1;
    double oy = x2;

    double best = INF;

    double bx = x1;
    double by = x2;

    // Проєкції на окремі стіни

    for (int i = 0; i < (int)walls.size(); i++)
    {
        double g = G(walls[i], ox, oy);

        double n2 =
            walls[i].a1 * walls[i].a1 +
            walls[i].a2 * walls[i].a2;

        double px =
            ox - g / n2 * walls[i].a1;

        double py =
            oy - g / n2 * walls[i].a2;

        if (feasible(px, py, walls, 1e-8))
        {
            double dist =
                (px - ox) * (px - ox) +
                (py - oy) * (py - oy);

            if (dist < best)
            {
                best = dist;

                bx = px;
                by = py;
            }
        }
    }

    // Перетини двох стін

    for (int i = 0; i < (int)walls.size(); i++)
    {
        for (int j = i + 1;
             j < (int)walls.size(); j++)
        {
            double det =
                walls[i].a1 * walls[j].a2 -
                walls[i].a2 * walls[j].a1;

            if (fabs(det) < 1e-12)
                continue;

            double px =
                (walls[i].b * walls[j].a2 -
                 walls[i].a2 * walls[j].b) / det;

            double py =
                (walls[i].a1 * walls[j].b -
                 walls[i].b * walls[j].a1) / det;

            if (feasible(px, py, walls, 1e-8))
            {
                double dist =
                    (px - ox) * (px - ox) +
                    (py - oy) * (py - oy);

                if (dist < best)
                {
                    best = dist;

                    bx = px;
                    by = py;
                }
            }
        }
    }

    x1 = bx;
    x2 = by;
}


// Активні стіни 

vector<int> activeSet(double x1, double x2,
                      const vector<Wall> &walls)
{
    vector<int> I;

    for (int i = 0;
         i < (int)walls.size(); i++)
    {
        if (fabs(G(walls[i],
                   x1, x2)) < ACTIVE_EPS)
        {
            I.push_back(i);
        }
    }

    return I;
}

bool inActive(const vector<int> &I,
              int j)
{
    for (int x : I)
    {
        if (x == j)
            return true;
    }

    return false;
}


// d = P grad F


void direction(double x1, double x2,
               const vector<Wall> &walls,
               const vector<int> &I,
               double &d1, double &d2)
{
    double g1, g2;

    gradF(x1, x2, g1, g2);

    // Немає активних стін
    if (I.empty())
    {
        d1 = -g1;
        d2 = -g2;

        return;
    }

    // Одна активна стіна
    if (I.size() == 1)
    {
        const Wall &a = walls[I[0]];

        double n2 =
            a.a1 * a.a1 +
            a.a2 * a.a2;

        double c =
            (a.a1 * g1 +
             a.a2 * g2) / n2;

        d1 = -(g1 - c * a.a1);
        d2 = -(g2 - c * a.a2);

        return;
    }

    // Дві незалежні стіни у R^2:
    // рухатись нікуди

    const Wall &a = walls[I[0]];
    const Wall &b = walls[I[1]];

    double det =
        a.a1 * b.a2 -
        a.a2 * b.a1;

    if (fabs(det) > 1e-12)
    {
        d1 = 0;
        d2 = 0;

        return;
    }
}


// Множники u


vector<double> multipliers(
    double x1, double x2,
    const vector<Wall> &walls,
    const vector<int> &I)
{
    vector<double> u;

    double g1, g2;

    gradF(x1, x2, g1, g2);

    if (I.empty())
        return u;

    // Одна стіна

    if (I.size() == 1)
    {
        const Wall &a = walls[I[0]];

        double n2 =
            a.a1 * a.a1 +
            a.a2 * a.a2;

        double value =
            -(a.a1 * g1 +
              a.a2 * g2) / n2;

        u.push_back(value);

        return u;
    }

    // Дві стіни

    const Wall &a = walls[I[0]];
    const Wall &b = walls[I[1]];

    double m11 =
        a.a1 * a.a1 +
        a.a2 * a.a2;

    double m12 =
        a.a1 * b.a1 +
        a.a2 * b.a2;

    double m22 =
        b.a1 * b.a1 +
        b.a2 * b.a2;

    double y1 =
        a.a1 * g1 +
        a.a2 * g2;

    double y2 =
        b.a1 * g1 +
        b.a2 * g2;

    double det =
        m11 * m22 -
        m12 * m12;

    double u1 =
        -(m22 * y1 -
          m12 * y2) / det;

    double u2 =
        -(-m12 * y1 +
          m11 * y2) / det;

    u.push_back(u1);
    u.push_back(u2);

    return u;
}


// F(x + h*d)


double lineF(double x1, double x2,
             double d1, double d2,
             double h)
{
    return F(x1 + h * d1,
             x2 + h * d2);
}


// Свенн + золотий переріз з парканом hmax


double lineSearch(double x1, double x2,
                  double d1, double d2,
                  double hmax)
{
    double delta = 0.01;

    double f0 =
        lineF(x1, x2,
              d1, d2, 0);

    double h1 = delta;

    if (isfinite(hmax) &&
        h1 > hmax)
    {
        h1 = hmax;
    }

    double f1 =
        lineF(x1, x2,
              d1, d2, h1);

    double a = 0;
    double b = h1;

    if (f1 < f0 &&
        (!isfinite(hmax) ||
         h1 < hmax - 1e-12))
    {
        double hPrev = 0;
        double h = h1;

        int k = 1;

        while (true)
        {
            double hNext =
                h + pow(2.0, k) * delta;

            if (isfinite(hmax) &&
                hNext > hmax)
            {
                hNext = hmax;
            }

            double fh =
                lineF(x1, x2,
                      d1, d2, h);

            double fn =
                lineF(x1, x2,
                      d1, d2, hNext);

            if (fn >= fh ||
                (isfinite(hmax) &&
                 fabs(hNext - hmax) <
                 1e-12))
            {
                a = hPrev;
                b = hNext;

                break;
            }

            hPrev = h;
            h = hNext;

            k++;

            if (k > 50)
            {
                a = hPrev;
                b = h;

                break;
            }
        }
    }

    // Золотий переріз

    double q =
        (sqrt(5.0) - 1) / 2;

    double eps = 1e-6;

    double left = a;
    double right = b;

    double xL =
        right -
        q * (right - left);

    double xR =
        left +
        q * (right - left);

    double fL =
        lineF(x1, x2,
              d1, d2, xL);

    double fR =
        lineF(x1, x2,
              d1, d2, xR);

    while (right - left > eps)
    {
        if (fL < fR)
        {
            right = xR;

            xR = xL;
            fR = fL;

            xL =
                right -
                q * (right - left);

            fL =
                lineF(x1, x2,
                      d1, d2, xL);
        }
        else
        {
            left = xL;

            xL = xR;
            fL = fR;

            xR =
                left +
                q * (right - left);

            fR =
                lineF(x1, x2,
                      d1, d2, xR);
        }
    }

    double h =
        (left + right) / 2;

    // Якщо мінімум саме на паркані

    if (isfinite(hmax))
    {
        double fh =
            lineF(x1, x2,
                  d1, d2, h);

        double fmax =
            lineF(x1, x2,
                  d1, d2, hmax);

        if (fmax <= fh + 1e-9)
            h = hmax;
    }

    return h;
}


// Метод Розена

Result rosen(string title,
             vector<Wall> walls,
             double x1, double x2)
{
    NF = 0;
    Ngrad = 0;

    cout << "\n==============================\n";
    cout << title << endl;
    cout << "==============================\n";

    cout << "Start: ("
         << x1 << "; "
         << x2 << ")\n";

    // Крок 0

    if (!feasible(x1, x2, walls))
    {
        projectPoint(x1, x2,
                     walls);

        cout << "Step 0 projection -> ("
             << x1 << "; "
             << x2 << ")"
             << ", F = "
             << Fraw(x1, x2)
             << endl;
    }

    vector<int> I =
        activeSet(x1, x2,
                  walls);

    int k = 0;

    cout << left
         << setw(4)  << "k"
         << setw(25) << "x"
         << setw(12) << "I"
         << setw(25) << "d"
         << setw(12) << "hmax"
         << setw(12) << "h"
         << setw(14) << "F"
         << "NF"
         << endl;

    while (k < 50)
    {
        double d1, d2;

        direction(x1, x2,
                  walls, I,
                  d1, d2);

        double norm =
            sqrt(d1 * d1 +
                 d2 * d2);

     
        // d = 0

        if (norm <= EPS)
        {
            vector<double> u =
                multipliers(
                    x1, x2,
                    walls, I);

            bool stop = true;

            int pos = -1;

            double minU = INF;

            for (int i = 0;
                 i < (int)u.size(); i++)
            {
                if (u[i] < -1e-9)
                {
                    stop = false;

                    if (u[i] < minU)
                    {
                        minU = u[i];
                        pos = i;
                    }
                }
            }

            if (stop)
            {
                cout << "STOP: d = 0, u = (";

                for (int i = 0;
                     i < (int)u.size(); i++)
                {
                    if (i)
                        cout << "; ";

                    cout << u[i];
                }

                cout << ")\n";

                return {
                    x1, x2,
                    I, u,
                    k,
                    NF,
                    Ngrad
                };
            }

            // Відпускаємо стіну
            // з найбільш від'ємним u

            cout << "Release wall "
                 << walls[I[pos]].name
                 << ", u = "
                 << u[pos]
                 << endl;

            I.erase(I.begin() + pos);

            continue;
        }

        
        // Паркан hmax

        double hmax =
            numeric_limits<double>
            ::infinity();

        int jHit = -1;

        for (int j = 0;
             j < (int)walls.size(); j++)
        {
            if (inActive(I, j))
                continue;

            double ad =
                walls[j].a1 * d1 +
                walls[j].a2 * d2;

            if (ad > 1e-12)
            {
                double hj =
                    (walls[j].b -
                     walls[j].a1 * x1 -
                     walls[j].a2 * x2)
                    / ad;

                if (hj < hmax)
                {
                    hmax = hj;
                    jHit = j;
                }
            }
        }

        double h =
            lineSearch(
                x1, x2,
                d1, d2,
                hmax);

        string active = "{";

        for (int i = 0;
             i < (int)I.size(); i++)
        {
            if (i)
                active += ",";

            active +=
                to_string(I[i] + 1);
        }

        active += "}";

        cout << left
             << setw(4) << k
             << setw(25)
             << ("(" +
                 to_string(x1) +
                 "; " +
                 to_string(x2) +
                 ")")
             << setw(12)
             << active
             << setw(25)
             << ("(" +
                 to_string(d1) +
                 "; " +
                 to_string(d2) +
                 ")");

        if (isfinite(hmax))
            cout << setw(12)
                 << hmax;
        else
            cout << setw(12)
                 << "inf";

        cout << setw(12)
             << h

             << setw(14)
             << Fraw(x1, x2)

             << NF
             << endl;

        // Нова точка

        x1 += h * d1;
        x2 += h * d2;

        // Якщо уперлись у паркан 
        // нова стіна активна

        if (jHit != -1 &&
            isfinite(hmax) &&
            fabs(h - hmax) <
            1e-7)
        {
            double gj =
                G(walls[jHit],
                  x1, x2);

            double n2 =
                walls[jHit].a1 *
                walls[jHit].a1 +
                walls[jHit].a2 *
                walls[jHit].a2;

            x1 -=
                gj / n2 *
                walls[jHit].a1;

            x2 -=
                gj / n2 *
                walls[jHit].a2;

            if (!inActive(I,
                          jHit))
            {
                I.push_back(jHit);
            }
        }

        k++;
    }

    return {
        x1, x2,
        I,
        multipliers(
            x1, x2,
            walls, I),
        k,
        NF,
        Ngrad
    };
}

// Коло: проєкція точки

void circleMethod(double x1,
                  double x2,
                  double R)
{
    NF = 0;
    Ngrad = 0;

    double h = 0.1;

    int k = 0;

    cout << "Start: ("
         << x1 << "; "
         << x2 << ")\n";

    while (k < 10000)
    {
        double g1, g2;

        gradF(x1, x2,
              g1, g2);

        double y1 =
            x1 - h * g1;

        double y2 =
            x2 - h * g2;

        double norm =
            sqrt(y1 * y1 +
                 y2 * y2);

        // якщо потрапили
        // всередину кола 
        // повертаємо на межу

        if (norm < R)
        {
            y1 = y1 * R / norm;
            y2 = y2 * R / norm;
        }

        double step =
            sqrt(
                (y1 - x1) *
                (y1 - x1) +

                (y2 - x2) *
                (y2 - x2)
            );

        x1 = y1;
        x2 = y2;

        k++;

        if (step <= EPS)
            break;
    }

    cout << "x = ("
         << x1 << "; "
         << x2 << ")\n";

    cout << "F = "
         << Fraw(x1, x2)
         << endl;

    cout << "norm = "
         << sqrt(x1 * x1 +
                 x2 * x2)
         << endl;

    cout << "iter = "
         << k
         << endl;

    cout << "Ngrad = "
         << Ngrad
         << endl;
}

// MAIN

int main()
{
    cout << fixed
         << setprecision(6);

    // Варіант 15

    Wall w1 =
    {
        1, -1, 3, "1"
    };

    Wall w2 =
    {
        0, 1, -2, "2"
    };

    Wall w3 =
    {
        0, -1, 4, "3"
    };


    // 3. Прогін 1

    Result run1 =
        rosen(
            "3. RUN 1 - WALL 1",
            {w1},
            9, 0
        );


    // 4. Прогін 2

    Result run2 =
        rosen(
            "4. RUN 2 - WALLS 1 + 2",
            {w1, w2},
            9, 0
        );


    // 4B. Демонстрація паркана hmax
    // дві стіни, допустимий старт

    Result fence =
        rosen(
            "4B. HMAX DEMO - WALLS 1 + 2",
            {w1, w2},
            -6, -7
        );


    // 5. Другий старт


    Result second =
        rosen(
            "5. SECOND START - WALL 1",
            {w1},
            -6, -7
        );


    // Демонстрація відпускання стіни


    Result release =
        rosen(
            "RELEASE WALL - WALLS 1 + 3",
            {w1, w3},
            -1, -4
        );


    // 6. Порівняння з ЛР 11


    cout << "\n==============================\n";
    cout << "6. COMPARISON WITH LAB 11\n";

    cout << left
         << setw(20) << "Method"
         << setw(10) << "steps"
         << setw(10) << "NF"
         << setw(10) << "Ngrad"
         << setw(26) << "x"
         << setw(14) << "F"
         << setw(14) << "porush"
         << "lambda"
         << endl;

    // ЛР 11, 1 обмеження

    cout << setw(20)
         << "Penalty LR11"

         << setw(10)
         << 5

         << setw(10)
         << 395

         << setw(10)
         << 18

         << setw(26)
         << "(1.666647;-1.333720)"

         << setw(14)
         << -34.002688

         << setw(14)
         << 0.000367

         << "7.332619"
         << endl;

    // Розен, 1 обмеження

    cout << setw(20)
         << "Rosen wall1"

         << setw(10)
         << run1.steps

         << setw(10)
         << run1.NF

         << setw(10)
         << run1.Ngrad;

    cout << setw(26)
         << ("(" +
             to_string(run1.x1) +
             ";" +
             to_string(run1.x2) +
             ")");

    cout << setw(14)
         << Fraw(run1.x1,
                 run1.x2);

    cout << setw(14)
         << max(
             0.0,
             G(w1,
               run1.x1,
               run1.x2));

    if (!run1.u.empty())
        cout << run1.u[0];

    cout << endl;

    // ЛР 11, два обмеження

    cout << setw(20)
         << "Penalty 2 LR11"

         << setw(10)
         << 5

         << setw(10)
         << 522

         << setw(10)
         << 23

         << setw(26)
         << "(1.001299;-1.999401)"

         << setw(14)
         << -30.016972

         << setw(14)
         << 0.000699

         << "(13.988;11.982)"
         << endl;

    // Розен, два обмеження

    cout << setw(20)
         << "Rosen walls1+2"

         << setw(10)
         << run2.steps

         << setw(10)
         << run2.NF

         << setw(10)
         << run2.Ngrad;

    cout << setw(26)
         << ("(" +
             to_string(run2.x1) +
             ";" +
             to_string(run2.x2) +
             ")");

    cout << setw(14)
         << Fraw(run2.x1,
                 run2.x2);

    double por =
        max(
            max(
                0.0,
                G(w1,
                  run2.x1,
                  run2.x2)
            ),
            max(
                0.0,
                G(w2,
                  run2.x1,
                  run2.x2)
            )
        );

    cout << setw(14)
         << por;

    cout << "(";

    for (int i = 0;
         i < (int)run2.u.size(); i++)
    {
        if (i)
            cout << ";";

        cout << run2.u[i];
    }

    cout << ")\n";


    // 7. Точки траєкторій для графіка


    cout << "\n==============================\n";
    cout << "7. TRAJECTORIES\n";


    cout << "Trajectory Rosen is printed "
         << "in RUN 1 and SECOND START.\n";

    cout << "LR11 trajectory:\n";

    cout << "(2.399999,-2.250000) -> "
         << "(1.811863,-1.514833) -> "
         << "(1.682685,-1.353540) -> "
         << "(1.668294,-1.335368) -> "
         << "(1.666647,-1.333720)\n";


    // 8. Коло, бонус


    cout << "\n==============================\n";
    cout << "8. CIRCLE R = 6\n";


    circleMethod(
        9, 0, 6
    );

    cout << endl;

    circleMethod(
        -9, 0, 6
    );

    return 0;
}