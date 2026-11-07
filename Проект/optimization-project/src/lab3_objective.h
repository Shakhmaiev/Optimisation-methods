
#pragma once

#include "image.h"
#include "metrics.h"
#include "filters.h"

#include <algorithm>
#include <functional>

using namespace std;

// Лічильник звернень до функції двох змінних
struct Oracle2 {

    function<double(double, double)> f;
    int n_f = 0;

    Oracle2(function<double(double, double)> func) {
        f = func;
    }

    double value(double a, double b) {
        n_f++;
        return f(a, b);
    }

    double peek(double a, double b) {
        return f(a, b);
    }

    void reset() {
        n_f = 0;
    }
};

// Цільова функція для білатерального фільтра
double bilateral_objective(
    const Image& noisy,
    const Image& clean,
    double sigma_s,
    double sigma_r
) {
    // Допустимі межі
    double s = clamp(sigma_s, 0.3, 5.0);
    double r = clamp(sigma_r, 0.01, 1.0);

    Image restored = bilateral_filter(noisy, s, r);

    double error = mse(restored, clean);

    // Штраф за вихід за межі
    double penalty = 1000.0 * (
        (sigma_s - s) * (sigma_s - s) +
        (sigma_r - r) * (sigma_r - r)
    );

    return error + penalty;
}
