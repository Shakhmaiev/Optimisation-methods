
#pragma once

#include "lab3_objective.h"
#include <array>
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>

using namespace std;

struct NMRow {
    int k;
    string operation;
    double a, b, f, size;
    int calls;
};

struct NMResult {
    double a, b, f;
    int iterations;
    int calls;
    vector<NMRow> rows;
};

inline NMResult nelder_mead(Oracle2& oracle,
                            double a0, double b0,
                            double step,
                            double tol = 1e-7) {

    auto clampA = [](double x) {
        return clamp(x, 0.3, 5.0);
    };

    auto clampB = [](double x) {
        return clamp(x, 0.01, 1.0);
    };

    // Три вершини початкового трикутника
    array<double, 3> a = {
        a0, clampA(a0 + step), a0
    };

    array<double, 3> b = {
        b0, b0, clampB(b0 + 0.4 * step)
    };

    array<double, 3> f;

    for (int i = 0; i < 3; i++)
        f[i] = oracle.value(a[i], b[i]);

    vector<NMRow> rows;
    int iterations = 0;

    for (int k = 1; k <= 200; k++) {

        // Сортуємо вершини від кращої до гіршої
        array<int, 3> order = {0, 1, 2};

        sort(order.begin(), order.end(), [&](int i, int j) {
            return f[i] < f[j];
        });

        auto oldA = a;
        auto oldB = b;
        auto oldF = f;

        for (int i = 0; i < 3; i++) {
            a[i] = oldA[order[i]];
            b[i] = oldB[order[i]];
            f[i] = oldF[order[i]];
        }

        // Центр двох кращих вершин
        double ca = (a[0] + a[1]) / 2.0;
        double cb = (b[0] + b[1]) / 2.0;

        // Відображення
        double ra = clampA(ca + (ca - a[2]));
        double rb = clampB(cb + (cb - b[2]));

        double fr = oracle.value(ra, rb);
        string operation;

        if (fr < f[0]) {

            // Розтягування
            double ea = clampA(ca + 2.0 * (ra - ca));
            double eb = clampB(cb + 2.0 * (rb - cb));

            double fe = oracle.value(ea, eb);

            if (fe < fr) {
                a[2] = ea;
                b[2] = eb;
                f[2] = fe;
                operation = "Expansion";
            } else {
                a[2] = ra;
                b[2] = rb;
                f[2] = fr;
                operation = "Reflection";
            }

        } else if (fr < f[1]) {

            a[2] = ra;
            b[2] = rb;
            f[2] = fr;
            operation = "Reflection";

        } else {

            // Стиснення
            double sa, sb;

            if (fr < f[2]) {
                sa = ca + 0.5 * (ra - ca);
                sb = cb + 0.5 * (rb - cb);
                operation = "Outside contraction";
            } else {
                sa = ca + 0.5 * (a[2] - ca);
                sb = cb + 0.5 * (b[2] - cb);
                operation = "Inside contraction";
            }

            sa = clampA(sa);
            sb = clampB(sb);

            double fs = oracle.value(sa, sb);

            if (fs < min(fr, f[2])) {
                a[2] = sa;
                b[2] = sb;
                f[2] = fs;
            } else {

                // Редукція трикутника
                for (int i = 1; i < 3; i++) {
                    a[i] = clampA(a[0] + 0.5 * (a[i] - a[0]));
                    b[i] = clampB(b[0] + 0.5 * (b[i] - b[0]));
                    f[i] = oracle.value(a[i], b[i]);
                }

                operation = "Reduction";
            }
        }

        // Знаходимо найкращу вершину
        int best = int(min_element(f.begin(), f.end()) - f.begin());

        double size = 0.0;

        for (int i = 0; i < 3; i++) {
            double da = a[i] - a[best];
            double db = b[i] - b[best];

            size = max(size, sqrt(da * da + db * db));
        }

        double spread = *max_element(f.begin(), f.end())
                      - *min_element(f.begin(), f.end());

        rows.push_back({
            k, operation,
            a[best], b[best], f[best],
            size, oracle.n_f
        });

        iterations = k;

        if (size <= tol && spread <= tol)
            break;
    }

    int best = int(min_element(f.begin(), f.end()) - f.begin());

    return {
        a[best], b[best], f[best],
        iterations, oracle.n_f, rows
    };
}
