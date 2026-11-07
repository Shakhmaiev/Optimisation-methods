
#pragma once

#include "operators.h"
#include <iostream>
#include <iomanip>
#include <random>
#include <cmath>
#include <algorithm>

inline void lab4_tests() {

    std::cout << "\nLAB 4 - OPERATORS TEST\n";

    std::mt19937 rng(12345);
    std::normal_distribution<double> normal(0.0, 1.0);
    std::uniform_real_distribution<double> uniform(0.0, 1.0);

    // TEST 1: Спряженість градієнта і дивергенції

    Image u(9, 7);
    Image px(9, 7);
    Image py(9, 7);

    for (size_t k = 0; k < u.px.size(); k++) {
        u.px[k] = normal(rng);
        px.px[k] = normal(rng);
        py.px[k] = normal(rng);
    }

    Field g = grad_xy(u);
    Image div = divergence(px, py);

    double left = 0.0;
    double right = 0.0;

    for (size_t k = 0; k < u.px.size(); k++) {

        left += g.x.px[k] * px.px[k]
              + g.y.px[k] * py.px[k];

        right -= u.px[k] * div.px[k];
    }

    double adjoint_error = std::abs(left - right);

    std::cout << "TEST 1 - ADJOINT CHECK\n";
    std::cout << "Left = " << left << "\n";
    std::cout << "Right = " << right << "\n";
    std::cout << "Difference = " << adjoint_error << "\n";

    // TEST 2: Аналітичний градієнт проти числового

    Image I(16, 16);
    Image noisy(16, 16);
    Image direction(16, 16);

    for (size_t k = 0; k < I.px.size(); k++) {
        I.px[k] = uniform(rng);
        noisy.px[k] = uniform(rng);
        direction.px[k] = normal(rng);
    }

    double lambda = 0.6548;
    double t = 1e-6;

    Image plus = I;
    Image minus = I;

    for (size_t k = 0; k < I.px.size(); k++) {

        plus.px[k] += t * direction.px[k];
        minus.px[k] -= t * direction.px[k];
    }

    double numeric = (
        energy(plus, noisy, lambda)
        - energy(minus, noisy, lambda)
    ) / (2.0 * t);

    Image gradient = energy_grad(I, noisy, lambda);

    double analytic = 0.0;

    for (size_t k = 0; k < I.px.size(); k++)
        analytic += gradient.px[k] * direction.px[k];

    double relative_error =
        std::abs(numeric - analytic)
        / std::max(std::abs(analytic), 1e-12);

    std::cout << "\nTEST 2 - GRADIENT CHECK\n";
    std::cout << "Numeric = " << numeric << "\n";
    std::cout << "Analytic = " << analytic << "\n";
    std::cout << "Relative error = " << relative_error << "\n";

    // Підсумок
    bool test1 = adjoint_error < 1e-10;
    bool test2 = relative_error < 1e-6;

    std::cout << "\nTEST RESULTS\n";
    std::cout << "Adjoint: " << (test1 ? "PASS" : "FAIL") << "\n";
    std::cout << "Gradient: " << (test2 ? "PASS" : "FAIL") << "\n";
}
