
#pragma once

#include "operators.h"
#include <cmath>

struct ReferenceResult {
    Image image;
    double E;
    double grad_norm;
    int iterations;
    bool converged;
};

inline ReferenceResult find_reference(
    const Image& noisy,
    double lambda,
    int max_iter = 3000,
    double tol = 1e-10
) {
    double L = 1.0 + 8.0 * lambda;
    double h = 1.0 / L;

    Image I = noisy;

    int iterations = 0;
    bool converged = false;

    for (int k = 0; k < max_iter; k++) {

        Image g = energy_grad(I, noisy, lambda);
        double g_norm = norm(g);

        if (g_norm < tol) {
            converged = true;
            break;
        }

        for (size_t i = 0; i < I.px.size(); i++) {
            I.px[i] -= h * g.px[i];
        }

        iterations++;
    }

    double final_grad = norm(
        energy_grad(I, noisy, lambda)
    );

    if (final_grad < tol)
        converged = true;

    return {
        I,
        energy(I, noisy, lambda),
        final_grad,
        iterations,
        converged
    };
}
