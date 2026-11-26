#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>
#include <string>
#include <limits>
#include <algorithm>

using namespace std;

const double EPS = 1e-9;

// Pivot

void pivot(vector<vector<double>>& T,
           vector<int>& basis,
           int row,
           int col)
{
    int R = T.size();
    int C = T[0].size();

    double p = T[row][col];

    for (int j = 0; j < C; j++)
        T[row][j] /= p;

    for (int i = 0; i < R; i++)
    {
        if (i == row)
            continue;

        double k = T[i][col];

        if (fabs(k) < EPS)
            continue;

        for (int j = 0; j < C; j++)
            T[i][j] -= k * T[row][j];
    }

    basis[row] = col;
}

// Вибір ведучого стовпця

int chooseColumn(const vector<vector<double>>& T,
                 string rule)
{
    int m = T.size() - 1;
    int n = T[0].size() - 1;

    int col = -1;

    if (rule == "bland")
    {
        for (int j = 0; j < n; j++)
        {
            if (T[m][j] < -EPS)
            {
                col = j;
                break;
            }
        }
    }
    else
    {
        double best = -EPS;

        for (int j = 0; j < n; j++)
        {
            if (T[m][j] < best)
            {
                best = T[m][j];
                col = j;
            }
        }
    }

    return col;
}

// Вибір ведучого рядка

int chooseRow(const vector<vector<double>>& T,
              int col)
{
    int m = T.size() - 1;
    int rhs = T[0].size() - 1;

    double best =
        numeric_limits<double>::infinity();

    int row = -1;

    for (int i = 0; i < m; i++)
    {
        if (T[i][col] > EPS)
        {
            double q =
                T[i][rhs] /
                T[i][col];

            if (q < best - EPS)
            {
                best = q;
                row = i;
            }
        }
    }

    return row;
}

// Друк таблиці

void printTable(const vector<vector<double>>& T,
                const vector<int>& basis,
                const vector<string>& names,
                int k)
{
    int m = T.size() - 1;
    int C = T[0].size();

    cout << "\nTable " << k << "\n";

    cout << left
         << setw(8) << "basis";

    for (int j = 0; j < C - 1; j++)
        cout << setw(10) << names[j];

    cout << setw(10) << "b"
         << endl;

    for (int i = 0; i < m; i++)
    {
        cout << left
             << setw(8)
             << names[basis[i]];

        for (int j = 0; j < C; j++)
            cout << setw(10)
                 << T[i][j];

        cout << endl;
    }

    cout << left
         << setw(8) << "Z";

    for (int j = 0; j < C; j++)
        cout << setw(10)
             << T[m][j];

    cout << endl;
}

// Звичайний simplex

void runSimplex(vector<vector<double>> T,
                vector<int> basis,
                vector<string> names,
                string rule)
{
    int k = 0;

    printTable(T,
               basis,
               names,
               k);

    while (true)
    {
        int col =
            chooseColumn(T,
                         rule);

        if (col == -1)
        {
            cout << "\nOPTIMUM\n";
            break;
        }

        int row =
            chooseRow(T,
                      col);

        if (row == -1)
        {
            cout << "\nZ UNBOUNDED\n";
            return;
        }

        cout << "\nEnter: "
             << names[col]
             << ", leave: "
             << names[basis[row]]
             << ", pivot = "
             << T[row][col]
             << endl;

        pivot(T,
              basis,
              row,
              col);

        k++;

        printTable(T,
                   basis,
                   names,
                   k);
    }

    int rhs =
        T[0].size() - 1;

    double x1 = 0;
    double x2 = 0;

    for (int i = 0;
         i < (int)basis.size();
         i++)
    {
        if (basis[i] == 0)
            x1 = T[i][rhs];

        if (basis[i] == 1)
            x2 = T[i][rhs];
    }

    cout << "\nResult:\n";

    cout << "x1 = "
         << x1 << endl;

    cout << "x2 = "
         << x2 << endl;

    cout << "Z = "
         << T.back()[rhs]
         << endl;

    // Тіньові ціни
    cout << "Shadow prices: ";

    for (int j = 2;
         j <= 4;
         j++)
    {
        cout << T.back()[j];

        if (j != 4)
            cout << "; ";
    }

    cout << endl;
}

// Двофазний метод для частини 2

void twoPhase()
{
    cout << "\n==============================\n";
    cout << "7. TWO-PHASE METHOD\n";
    cout << "==============================\n";

    /*
    x1 + 2x2 >= 10
    x1 <= 8
    x2 <= 7

    x1 + 2x2 - s1 + a1 = 10
    x1 + s2 = 8
    x2 + s3 = 7
    */

    vector<string> names =
    {
        "x1", "x2",
        "s1", "s2", "s3",
        "a1"
    };

    vector<vector<double>> T =
    {
        { 1, 2, -1, 0, 0, 1, 10 },
        { 1, 0,  0, 1, 0, 0,  8 },
        { 0, 1,  0, 0, 1, 0,  7 },

        // Фаза 1:
        // max(-a1)
        // після приведення по базису a1
        { -1, -2, 1, 0, 0, 0, -10 }
    };

    vector<int> basis =
    {
        5, 3, 4
    };

    cout << "\nPHASE 1\n";

    int k = 0;

    printTable(T,
               basis,
               names,
               k);

    while (true)
    {
        int col =
            chooseColumn(T,
                         "dantzig");

        if (col == -1)
            break;

        int row =
            chooseRow(T,
                      col);

        if (row == -1)
        {
            cout << "No feasible solution\n";
            return;
        }

        cout << "\nEnter: "
             << names[col]
             << ", leave: "
             << names[basis[row]]
             << endl;

        pivot(T,
              basis,
              row,
              col);

        k++;

        printTable(T,
                   basis,
                   names,
                   k);
    }

    cout << "\nPhase 1 result = "
         << T.back().back()
         << endl;

    // Видаляємо artificial a1

    for (auto& row : T)
        row.erase(row.begin() + 5);

    names.erase(names.begin() + 5);

    // Фаза 2
    // max(-W) = -2x1 - 8x2

    T.back() =
    {
        2, 8, 0, 0, 0, 0
    };

    // приводимо рядок Z до поточного базису

    for (int i = 0;
         i < 3;
         i++)
    {
        int bj = basis[i];

        double q =
            T.back()[bj];

        if (fabs(q) > EPS)
        {
            for (int j = 0;
                 j < (int)T[0].size();
                 j++)
            {
                T.back()[j] -=
                    q * T[i][j];
            }
        }
    }

    cout << "\nPHASE 2\n";

    k = 0;

    printTable(T,
               basis,
               names,
               k);

    while (true)
    {
        int col =
            chooseColumn(T,
                         "dantzig");

        if (col == -1)
            break;

        int row =
            chooseRow(T,
                      col);

        if (row == -1)
        {
            cout << "Unbounded\n";
            return;
        }

        cout << "\nEnter: "
             << names[col]
             << ", leave: "
             << names[basis[row]]
             << endl;

        pivot(T,
              basis,
              row,
              col);

        k++;

        printTable(T,
                   basis,
                   names,
                   k);
    }

    int rhs =
        T[0].size() - 1;

    double x1 = 0;
    double x2 = 0;

    for (int i = 0;
         i < 3;
         i++)
    {
        if (basis[i] == 0)
            x1 = T[i][rhs];

        if (basis[i] == 1)
            x2 = T[i][rhs];
    }

    double maxMinusW =
        T.back()[rhs];

    double W =
        -maxMinusW;

    cout << "\nResult part 2:\n";

    cout << "x1 = "
         << x1 << endl;

    cout << "x2 = "
         << x2 << endl;

    cout << "W = "
         << W << endl;
}


// ПУНКТ 8. МЕТОД РОЗЕНА ДЛЯ ЛІНІЙНОЇ ЗАДАЧІ


struct RWall
{
    double a1, a2, b;
    string name;
};

double RG(const RWall& w,
          double x1, double x2)
{
    return w.a1 * x1 +
           w.a2 * x2 -
           w.b;
}

// Мінімізуємо F = -Z
void Rgrad(double& g1,
           double& g2)
{
    g1 = -2.0;
    g2 = -4.0;
}

bool RisActive(const vector<int>& I,
               int j)
{
    for (int q : I)
        if (q == j)
            return true;

    return false;
}

// d = -P * gradF
void Rdirection(double x1,
                double x2,
                const vector<RWall>& walls,
                const vector<int>& I,
                double& d1,
                double& d2)
{
    double g1, g2;

    Rgrad(g1, g2);

    if (I.empty())
    {
        d1 = -g1;
        d2 = -g2;
        return;
    }

    if (I.size() == 1)
    {
        const RWall& a =
            walls[I[0]];

        double n2 =
            a.a1 * a.a1 +
            a.a2 * a.a2;

        double q =
            (a.a1 * g1 +
             a.a2 * g2) / n2;

        d1 =
            -(g1 - q * a.a1);

        d2 =
            -(g2 - q * a.a2);

        return;
    }

    // У R2 дві незалежні активні стіни
    // повністю блокують рух

    const RWall& a =
        walls[I[0]];

    const RWall& b =
        walls[I[1]];

    double det =
        a.a1 * b.a2 -
        a.a2 * b.a1;

    if (fabs(det) > EPS)
    {
        d1 = 0;
        d2 = 0;
    }
}

// Множники активних стін
vector<double> Rmultipliers(
    const vector<RWall>& walls,
    const vector<int>& I)
{
    vector<double> u;

    double g1, g2;

    Rgrad(g1, g2);

    if (I.empty())
        return u;

    if (I.size() == 1)
    {
        const RWall& a =
            walls[I[0]];

        double n2 =
            a.a1 * a.a1 +
            a.a2 * a.a2;

        double v =
            -(a.a1 * g1 +
              a.a2 * g2) / n2;

        u.push_back(v);

        return u;
    }

    const RWall& a =
        walls[I[0]];

    const RWall& b =
        walls[I[1]];

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

void rosenLP()
{
    cout << "\n==============================\n";
    cout << "8. ROSEN METHOD\n";
    cout << "==============================\n";

    // Стіни у формі a*x <= b

    vector<RWall> walls =
    {
        // x1 <= 3
        { 1,  0,  3, "x1 <= 3" },

        // x2 <= 8
        { 0,  1,  8, "x2 <= 8" },

        // 4x1 + 3x2 <= 28
        { 4,  3, 28, "4x1+3x2 <= 28" },

        // x1 >= 0 -> -x1 <= 0
        {-1,  0,  0, "x1 >= 0" },

        // x2 >= 0 -> -x2 <= 0
        { 0, -1,  0, "x2 >= 0" }
    };

    double x1 = 0;
    double x2 = 0;

    // У старті активні:
    // x1 = 0 та x2 = 0

    vector<int> I =
    {
        3, 4
    };

    int k = 0;

    cout << fixed
         << setprecision(6);

    cout << left
         << setw(5)  << "k"
         << setw(24) << "x"
         << setw(24) << "active"
         << setw(24) << "d"
         << setw(12) << "hmax"
         << setw(12) << "Z"
         << endl;

    while (k < 20)
    {
        double d1, d2;

        Rdirection(
            x1, x2,
            walls, I,
            d1, d2
        );

        double norm =
            sqrt(d1 * d1 +
                 d2 * d2);


        // Якщо d = 0 -> дивимось множники

        if (norm < 1e-9)
        {
            vector<double> u =
                Rmultipliers(
                    walls, I);

            bool stop = true;

            int pos = -1;
            double minU = 1e100;

            for (int i = 0;
                 i < (int)u.size();
                 i++)
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
                cout << "\nSTOP\n";

                cout << "x* = ("
                     << x1 << "; "
                     << x2 << ")\n";

                cout << "Z* = "
                     << 2 * x1 +
                        4 * x2
                     << endl;

                cout << "u = (";

                for (int i = 0;
                     i < (int)u.size();
                     i++)
                {
                    if (i)
                        cout << "; ";

                    cout << u[i];
                }

                cout << ")\n";

                cout << "Active walls:\n";

                for (int j : I)
                    cout << "  "
                         << walls[j].name
                         << endl;

                break;
            }

            cout << "Release: "
                 << walls[I[pos]].name
                 << ", u = "
                 << u[pos]
                 << endl;

            I.erase(
                I.begin() + pos
            );

            continue;
        }


        // hmax


        double hmax =
            numeric_limits<double>
            ::infinity();

        int jhit = -1;

        for (int j = 0;
             j < (int)walls.size();
             j++)
        {
            if (RisActive(I, j))
                continue;

            double ad =
                walls[j].a1 * d1 +
                walls[j].a2 * d2;

            if (ad > 1e-12)
            {
                double h =
                    (walls[j].b -
                     walls[j].a1 * x1 -
                     walls[j].a2 * x2)
                    / ad;

                if (h < hmax)
                {
                    hmax = h;
                    jhit = j;
                }
            }
        }

        // Для лінійної цілі вздовж d
        // Z постійно росте.
        // Тому завжди йдемо до паркана.

        double h = hmax;

        string act = "{";

        for (int i = 0;
             i < (int)I.size();
             i++)
        {
            if (i)
                act += ",";

            act +=
                to_string(I[i] + 1);
        }

        act += "}";

        cout << left
             << setw(5) << k

             << setw(24)
             << ("(" +
                 to_string(x1) +
                 ";" +
                 to_string(x2) +
                 ")")

             << setw(24)
             << act

             << setw(24)
             << ("(" +
                 to_string(d1) +
                 ";" +
                 to_string(d2) +
                 ")")

             << setw(12)
             << hmax

             << setw(12)
             << 2 * x1 +
                4 * x2

             << endl;

        x1 += h * d1;
        x2 += h * d2;

        // Явно ставимо точку на нову стіну

        if (jhit != -1)
        {
            double gj =
                RG(walls[jhit],
                   x1, x2);

            double n2 =
                walls[jhit].a1 *
                walls[jhit].a1 +

                walls[jhit].a2 *
                walls[jhit].a2;

            x1 -=
                gj / n2 *
                walls[jhit].a1;

            x2 -=
                gj / n2 *
                walls[jhit].a2;

            I.push_back(jhit);
        }

        k++;
    }

    cout << "\nSimplex Dantzig path:\n";
    cout << "(0;0) -> (0;8) -> (1;8)\n";

    cout << "Rosen result must be the same.\n";

    cout << "Shadow prices from simplex: "
         << "(0; 2.5; 0.5)\n";
}

// MAIN

int main()
{
    cout << fixed
         << setprecision(6);

    // 1. VARIANT 15

    cout << "==============================\n";
    cout << "1. VARIANT 15\n";
    cout << "==============================\n";

    cout << "Z = 2*x1 + 4*x2 -> max\n";
    cout << "x1 <= 3\n";
    cout << "x2 <= 8\n";
    cout << "4*x1 + 3*x2 <= 28\n";

    // 2. VERTICES

    cout << "\n==============================\n";
    cout << "2. VERTICES\n";
    cout << "==============================\n";

    cout << "(0;0)     Z = 0\n";
    cout << "(3;0)     Z = 6\n";
    cout << "(3;5.333333) Z = 27.333333\n";
    cout << "(1;8)     Z = 34\n";
    cout << "(0;8)     Z = 32\n";

    // 3-4. DANTZIG

    vector<string> names =
    {
        "x1", "x2",
        "s1", "s2", "s3"
    };

    vector<vector<double>> T =
    {
        { 1, 0, 1, 0, 0, 3 },
        { 0, 1, 0, 1, 0, 8 },
        { 4, 3, 0, 0, 1, 28 },

        { -2, -4, 0, 0, 0, 0 }
    };

    vector<int> basis =
    {
        2, 3, 4
    };

    cout << "\n==============================\n";
    cout << "3-4. DANTZIG SIMPLEX\n";
    cout << "==============================\n";

    runSimplex(T,
               basis,
               names,
               "dantzig");

    // 5. BLAND

    cout << "\n==============================\n";
    cout << "5. BLAND RULE\n";
    cout << "==============================\n";

    runSimplex(T,
               basis,
               names,
               "bland");

    // 6. DUAL CHECK

    cout << "\n==============================\n";
    cout << "6. SHADOW PRICES / DUAL\n";
    cout << "==============================\n";

    cout << "y = (0; 2.5; 0.5)\n";

    cout << "Dual:\n";

    cout << "W = 3*y1 + 8*y2 + 28*y3 -> min\n";

    cout << "y1 + 4*y3 >= 2\n";

    cout << "y2 + 3*y3 >= 4\n";

    cout << "y >= 0\n";

    cout << "Check:\n";

    cout << "W = "
         << 3 * 0
          + 8 * 2.5
          + 28 * 0.5
         << endl;

    cout << "Z = 34\n";


    // 7. TWO PHASE

    twoPhase();


    // 8. ROSEN - BONUS

    rosenLP();


    // SUMMARY

    cout << "\n==============================\n";
    cout << "SUMMARY\n";
    cout << "==============================\n";

    cout << "Part 1:\n";
    cout << "x* = (1; 8)\n";
    cout << "Z* = 34\n";
    cout << "shadow prices = (0; 2.5; 0.5)\n";

    cout << "\nPart 2:\n";
    cout << "x* = (8; 1)\n";
    cout << "W* = 24\n";

    return 0;
}