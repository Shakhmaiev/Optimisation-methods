
#pragma once

#include "operators.h"
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>

struct GradRow {
    int k;
    double E;
    double grad_norm;
    double step;
    int calls;
    int grad_calls;
};

struct GradResult {
    Image image;
    int iterations;
    int calls;
    int grad_calls;
    double E;
    double grad_norm;
    bool converged;
    std::vector<GradRow> rows;
};

inline GradResult gradient_method(
    const Image& noisy,
    double lambda,
    const std::string& kind,
    double beta = 0.0,
    int max_iter = 200,
    double tol = 1e-6,
double custom_step = -1.0
) {
    
double L = 1.0 + 8.0 * lambda;
double mu = 1.0;

double h = 1.0 / L;

// Оптимальний крок методу важкої кульки
// Use fixed baseline step 1/L for Heavy Ball, as in the practical exercise.

if (custom_step > 0.0) {
    h = custom_step;
}
    Image I = noisy;
    Image previous = noisy;

    int calls = 0;
    int grad_calls = 0;
    int iterations = 0;
    bool converged = false;

    std::vector<GradRow> rows;

    double trial = h;

    for (int k = 0; k < max_iter; k++) {

        Image Y = I;

        // Інерційна точка для Нестерова
        if (kind == "nag") {
            double b = beta;

            if (beta < 0.0)
                b = double(k) / double(k + 3);

            for (size_t i = 0; i < I.px.size(); i++) {
                Y.px[i] = I.px[i]
                        + b * (I.px[i] - previous.px[i]);
            }
        }

        Image g = energy_grad(Y, noisy, lambda);
        grad_calls++;

        double gn = norm(g);

        double E = energy(I, noisy, lambda);
        calls++;

        rows.push_back({k, E, gn, h, calls, grad_calls});

        if (gn <= tol) {
            converged = true;
            break;
        }

        Image next = I;

        if (kind == "armijo") {

            double g2 = gn * gn;
            double s = trial;

            for (int t = 0; t < 40; t++) {

                for (size_t i = 0; i < I.px.size(); i++)
                    next.px[i] = I.px[i] - s * g.px[i];

                double E_new = energy(next, noisy, lambda);
                calls++;

                if (E_new <= E - 1e-4 * s * g2)
                    break;

                s *= 0.5;
            }

            h = s;
            trial = std::min(s * 1.5, 2.0 / L);

        } else if (kind == "hb") {

            // Метод важкої кульки
            for (size_t i = 0; i < I.px.size(); i++) {

                next.px[i] = I.px[i]
                    - h * g.px[i]
                    + beta * (I.px[i] - previous.px[i]);
            }

        } else if (kind == "nag") {

            // Метод Нестерова
            for (size_t i = 0; i < I.px.size(); i++)
                next.px[i] = Y.px[i] - h * g.px[i];

        } else {

            // Звичайний градієнтний спуск
            for (size_t i = 0; i < I.px.size(); i++)
                next.px[i] = I.px[i] - h * g.px[i];
        }

        previous = I;
        I = next;
        iterations++;
    }

    double final_grad = norm(energy_grad(I, noisy, lambda));
    grad_calls++;

    double final_E = energy(I, noisy, lambda);
    calls++;

    if (final_grad <= tol)
        converged = true;

    return {
        I, iterations, calls, grad_calls,
        final_E, final_grad, converged, rows
    };
}
