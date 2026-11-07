
#pragma once

#include "oracle.h"
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>

using namespace std;

// Один рядок таблиці
struct SvennRow {
    int k;
    double x;
    double f;
    double delta;
    int p;
    string comment;
};

// Результати методу
struct SvennResult {
    double a, b;
    int iterations;
    int p;
    int n_f;
    vector<SvennRow> rows;
};

SvennResult svenn(Oracle& oracle, double x0, double delta,
                  double LO = 0.05, double HI = 8.0) {


    vector<SvennRow> rows;

    // Обчислення функції та запис у таблицю
    auto ask = [&](double x, int k, double step,
                   int p, string comment) {

        double value = oracle.value(x);

        rows.push_back({k, x, value, step, p, comment});

        return value;
    };

    // Початкові точки
    double xl = clamp(x0 - delta, LO, HI);
    double xr = clamp(x0 + delta, LO, HI);

    double f0 = ask(x0, 0, delta, 0, "Start");
    double fr = ask(xr, 0, delta, 0, "Right");
    double fl = ask(xl, 0, -delta, 1, "Left");

    // Мінімум між пробними точками
    if (fl >= f0 && fr >= f0) {
        return {xl, xr, 0, 0, oracle.n_f, rows};
    }

    // Визначення напрямку
    double step;
    int p;

    if (fr <= f0) {
        step = delta;
        p = 0;
    } else {
        step = -delta;
        p = 1;
    }

    double previous = x0;
    double current = (p == 0) ? xr : xl;
    double fcurrent = (p == 0) ? fr : fl;

    // Подвоєння кроку
    for (int k = 1; k <= 40; k++) {

        double dk = pow(2.0, k) * step;
        double next = clamp(current + dk, LO, HI);

        // Досягнення межі
        if (next == current) {
            return {min(previous, current),
                    max(previous, current),
                    k, p, oracle.n_f, rows};
        }

        double fnext = ask(next, k, dk, p, "Expansion");

        // Функція почала зростати
        if (fnext >= fcurrent) {

            rows.back().comment = "Stop";

            return {min(previous, next),
                    max(previous, next),
                    k, p, oracle.n_f, rows};
        }

        previous = current;
        current = next;
        fcurrent = fnext;
    }

    return {min(previous, current),
            max(previous, current),
            40, p, oracle.n_f, rows};
}
