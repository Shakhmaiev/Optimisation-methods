
#pragma once

#include "oracle.h"
#include <cmath>
#include <vector>

using namespace std;

// Один рядок таблиці
struct SearchRow {
    int k;
    double a, b;
    double x1, x2;
    double f1, f2;
    double length;
};

// Результат методу
struct SearchResult {
    double x;
    double f;
    double a, b;
    int iterations;
    int n_f;
    vector<SearchRow> rows;
};

// Метод дихотомії
SearchResult dichotomy(Oracle& oracle,
                       double a, double b, double eps) {

    double delta = eps / 3.0;
    int k = 0;

    vector<SearchRow> rows;

    while (b - a > eps && k < 500) {

        double x1 = (a + b - delta) / 2.0;
        double x2 = (a + b + delta) / 2.0;

        double f1 = oracle.value(x1);
        double f2 = oracle.value(x2);

        if (f1 <= f2) {
            b = x2;
        } else {
            a = x1;
        }

        k++;

        rows.push_back({
            k, a, b, x1, x2, f1, f2, b - a
        });
    }

    double x = (a + b) / 2.0;
    double f = oracle.value(x);

    return {x, f, a, b, k, oracle.n_f, rows};
}

// Метод золотого перерізу
SearchResult golden_section(Oracle& oracle,
                            double a, double b, double eps) {

    const double tau = (sqrt(5.0) - 1.0) / 2.0;

    double x1 = b - tau * (b - a);
    double x2 = a + tau * (b - a);

    double f1 = oracle.value(x1);
    double f2 = oracle.value(x2);

    int k = 0;
    vector<SearchRow> rows;

    while (b - a > eps && k < 500) {

        if (f1 <= f2) {
            b = x2;

            x2 = x1;
            f2 = f1;

            x1 = b - tau * (b - a);
            f1 = oracle.value(x1);
        } else {
            a = x1;

            x1 = x2;
            f1 = f2;

            x2 = a + tau * (b - a);
            f2 = oracle.value(x2);
        }

        k++;

        rows.push_back({
            k, a, b, x1, x2, f1, f2, b - a
        });
    }

    double x = (a + b) / 2.0;
    double f = oracle.value(x);

    return {x, f, a, b, k, oracle.n_f, rows};
}
