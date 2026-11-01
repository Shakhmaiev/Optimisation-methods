
#pragma once

#include "lab3_objective.h"
#include <algorithm>
#include <cmath>
#include <vector>

using namespace std;

struct PowellRow {
    int cycle;
    double a, b, f, shift;
    int calls;
};

struct PowellResult {
    double a, b, f;
    int cycles, calls;
    vector<PowellRow> rows;
};


struct LineSearchStats {
    int svenn = 0;
    int golden = 0;
    int final_check = 0;
};


// Пошук мінімуму вздовж напрямку
inline double line_search(Oracle2& oracle,
                          double& a, double& b,
                          double da, double db,
                          double eps = 1e-4,
LineSearchStats* stats = nullptr) {

    double start_a = a;
    double start_b = b;

    auto f = [&](double alpha) {
        double x = clamp(start_a + alpha * da, 0.3, 5.0);
        double y = clamp(start_b + alpha * db, 0.01, 1.0);
        return oracle.value(x, y);
    };

    int calls_start = oracle.n_f;
    // Метод Свенна
    double delta = 0.1;
    double f0 = f(0.0);
    double fr = f(delta);
    double fl = f(-delta);

    double left, right;

    if (f0 <= fr && f0 <= fl) {
        left = -delta;
        right = delta;
    } else {

        double step = (fr <= f0) ? delta : -delta;
        double previous = 0.0;
        double current = step;
        double fcurrent = (step > 0) ? fr : fl;

        left = -1000.0;
        right = 1000.0;

        for (int k = 1; k <= 40; k++) {

            double next = clamp(
                current + pow(2.0, k) * step,
                -1000.0, 1000.0
            );

            if (next == current) {
                left = min(previous, current);
                right = max(previous, current);
                break;
            }

            double fnext = f(next);

            if (fnext >= fcurrent) {
                left = min(previous, next);
                right = max(previous, next);
                break;
            }

            previous = current;
            current = next;
            fcurrent = fnext;
        }
    }

    int calls_after_svenn = oracle.n_f;
    // Золотий переріз
    const double tau = (sqrt(5.0) - 1.0) / 2.0;

    double x1 = right - tau * (right - left);
    double x2 = left + tau * (right - left);

    double f1 = f(x1);
    double f2 = f(x2);

    for (int k = 0; k < 200 && right - left > eps; k++) {

        if (f1 <= f2) {
            right = x2;
            x2 = x1;
            f2 = f1;
            x1 = right - tau * (right - left);
            f1 = f(x1);
        } else {
            left = x1;
            x1 = x2;
            f1 = f2;
            x2 = left + tau * (right - left);
            f2 = f(x2);
        }
    }
    double alpha = (left + right) / 2.0;

// Фінальне обчислення золотого перерізу
double best_line_value = f(alpha);

int calls_after_golden = oracle.n_f;

    a = clamp(start_a + alpha * da, 0.3, 5.0);
    b = clamp(start_b + alpha * db, 0.01, 1.0);

        double result = oracle.value(a, b);

    if (stats != nullptr) {
        stats->svenn = calls_after_svenn - calls_start;
        stats->golden = calls_after_golden - calls_after_svenn;
        stats->final_check = oracle.n_f - calls_after_golden;
    }

    return result;
}

// Метод Пауелла
inline PowellResult powell(Oracle2& oracle,
                            double a0, double b0,
                            double tol = 1e-4) {

    double a = a0;
    double b = b0;
    double fx = oracle.value(a, b);

    // Два початкові напрямки
    double da[2] = {1.0, 0.0};
    double db[2] = {0.0, 1.0};

    vector<PowellRow> rows;
    int cycles = 0;

    for (int k = 1; k <= 50; k++) {

        int cycle_calls_start = oracle.n_f;
        double old_a = a;
        double old_b = b;
        double old_f = fx;

        double gains[2];
        LineSearchStats first_cycle_stats[3];

        // Два координатні пошуки. Точка циклу фіксується до
        // додаткового пошуку за сумарним напрямком (як у методичці).
        for (int i = 0; i < 2; i++) {
            double previous_f = fx;
            fx = line_search(oracle, a, b, da[i], db[i], 1e-4,
                             &first_cycle_stats[i]);
            gains[i] = previous_f - fx;
        }

        double move_a = a - old_a;
        double move_b = b - old_b;
        double shift = hypot(move_a, move_b);
        cycles = k;
        if (shift < tol ||
            abs(old_f - fx) <= 1e-12 * (1.0 + abs(fx))) {
            rows.push_back({k, a, b, fx, shift, oracle.n_f});
            cout << "Powell cycle " << k << ": search1="
                 << first_cycle_stats[0].svenn + first_cycle_stats[0].golden + first_cycle_stats[0].final_check
                 << ", search2="
                 << first_cycle_stats[1].svenn + first_cycle_stats[1].golden + first_cycle_stats[1].final_check
                 << ", search3=SKIPPED, calls=" << oracle.n_f-cycle_calls_start
                 << ", reason=convergence" << endl;
            break;
        }

        // Нова координата напрямку
        double new_da = move_a / shift;
        double new_db = move_b / shift;

        fx = line_search(
    oracle, a, b, new_da, new_db,
    1e-4, &first_cycle_stats[2]
);

        if (k == 1) {
            cout << "\nLAB 3 - EXPERIMENT 2\n";

            int total = 0;

            for (int i = 0; i < 3; i++) {
                auto st = first_cycle_stats[i];

                int sum = st.svenn
                        + st.golden
                        + st.final_check;

                total += sum;

                cout << "Line search " << i + 1
                     << ": Svenn = " << st.svenn
                     << ", Golden = " << st.golden
                     << ", Final check = " << st.final_check
                     << ", Total = " << sum << endl;
            }

            cout << "Initial point evaluation = 1\n";
            cout << "First cycle total calls = "
                 << total + 1 << endl;
        }


        rows.push_back({k, a, b, fx, hypot(a-old_a, b-old_b), oracle.n_f});
        cout << "Powell cycle " << k << ": search1="
             << first_cycle_stats[0].svenn + first_cycle_stats[0].golden + first_cycle_stats[0].final_check
             << ", search2="
             << first_cycle_stats[1].svenn + first_cycle_stats[1].golden + first_cycle_stats[1].final_check
             << ", search3="
             << first_cycle_stats[2].svenn + first_cycle_stats[2].golden + first_cycle_stats[2].final_check
             << ", calls=" << oracle.n_f-cycle_calls_start << endl;

        // Замінюємо напрямок з найбільшим виграшем
        int index = (gains[0] >= gains[1]) ? 0 : 1;

        da[index] = new_da;
        db[index] = new_db;

        // Перезапуск напрямків кожні три цикли.

        if (k % 3 == 0) {
            da[0] = 1.0;
            db[0] = 0.0;

            da[1] = 0.0;
            db[1] = 1.0;
        }
    }

    return {a, b, fx, cycles, oracle.n_f, rows};
}
