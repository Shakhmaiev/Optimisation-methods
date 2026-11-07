
#pragma once

#include "operators.h"
#include "metrics.h"
#include <algorithm>

// Знаходимо мінімум енергії Тихонова
// градієнтним спуском зі сталим кроком
inline Image tikhonov_solve(
    const Image& noisy,
    double lambda,
    int iterations = 300
) {
    double L = 1.0 + 8.0 * lambda;
    double h = 1.0 / L;

    Image I = noisy;

    for (int k = 0; k < iterations; k++) {

        Image grad = energy_grad(I, noisy, lambda);

        for (size_t i = 0; i < I.px.size(); i++) {
            I.px[i] -= h * grad.px[i];
        }
    }

    return I;
}

// Зовнішня цільова функція для підбору lambda
inline double lambda_objective(
    const Image& noisy,
    const Image& clean,
    double lambda
) {
    double lam = std::clamp(lambda, 0.0001, 50.0);

    Image result = tikhonov_solve(noisy, lam);

    double penalty = 1000.0 *
        (lambda - lam) * (lambda - lam);

    return mse(result, clean) + penalty;
}
