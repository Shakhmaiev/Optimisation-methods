
#pragma once

#include "operators.h"
#include <algorithm>
#include <cmath>
#include <vector>

struct Violation {
    double max_value;
    int count;
};

inline Image project(
    Image I,
    double lower = 0.0,
    double upper = 1.0
) {
    for (double& value : I.px) {
        value = std::clamp(value, lower, upper);
    }

    return I;
}

inline Violation violation(
    const Image& I,
    double lower = 0.0,
    double upper = 1.0
) {
    Violation result{0.0, 0};

    for (double value : I.px) {

        double error = std::max(
            std::max(lower - value, value - upper),
            0.0
        );

        if (error > 0.0) {
            result.count++;
            result.max_value =
                std::max(result.max_value, error);
        }
    }

    return result;
}

inline int active_constraints(
    const Image& I,
    double eps = 1e-12
) {
    int count = 0;

    for (double value : I.px) {
        if (value <= eps || value >= 1.0 - eps)
            count++;
    }

    return count;
}


struct ProjectedRow {
    int k;
    double E;
    double G;
    int active;
    double violation;
    int calls;
    int grad_calls;
};

struct ProjectedResult {
    Image image;
    int iterations;
    int calls;
    int grad_calls;
    double E;
    double G;
    bool converged;
    std::vector<ProjectedRow> rows;
};

inline ProjectedResult projected_gradient(
    const Image& raw,
    double lambda,
    int max_iter = 200,
    double tol = 1e-6
) {
    double L = 1.0 + 8.0 * lambda;
    double h = 1.0 / L;

    Image I = project(raw);

    int calls = 0;
    int grad_calls = 0;
    int iterations = 0;

    double G = 0.0;
    bool converged = false;

    std::vector<ProjectedRow> rows;

    for (int k = 0; k <= max_iter; k++) {

        Image grad = energy_grad(I, raw, lambda);
        grad_calls++;

        Image next = I;

        for (size_t i = 0; i < I.px.size(); i++) {
            next.px[i] = std::clamp(
                I.px[i] - h * grad.px[i],
                0.0,
                1.0
            );
        }

        double sum = 0.0;

        for (size_t i = 0; i < I.px.size(); i++) {
            double d = (I.px[i] - next.px[i]) / h;
            sum += d * d;
        }

        G = std::sqrt(sum);

        double E = energy(I, raw, lambda);
        calls++;

        Violation v = violation(I);

        rows.push_back({
            k,
            E,
            G,
            active_constraints(I),
            v.max_value,
            calls,
            grad_calls
        });

        if (G <= tol) {
            converged = true;
            break;
        }

        if (k == max_iter)
            break;

        I = next;
        iterations++;
    }

    return {
        I,
        iterations,
        calls,
        grad_calls,
        rows.back().E,
        G,
        converged,
        rows
    };
}



struct PenaltyRow {
    int outer;
    double r;
    double Lr;
    double h;
    double E;
    double P;
    double rP;
    double max_violation;
    int outside;
    int inner_iterations;
    long calls;
    long grad_calls;
};

struct PenaltyResult {
    Image image;
    double E;
    double P;
    int outer_iterations;
    long calls;
    long grad_calls;
    std::vector<PenaltyRow> rows;
};

// Сума квадратів порушень
inline double penalty_value(const Image& I) {
    double sum = 0.0;

    for (double v : I.px) {
        double below = std::max(0.0, -v);
        double above = std::max(0.0, v - 1.0);

        sum += below * below + above * above;
    }

    return sum;
}

// Градієнт штрафної функції
inline Image penalty_grad(const Image& I) {
    Image g(I.h, I.w);

    for (size_t i = 0; i < I.px.size(); i++) {
        double v = I.px[i];

        if (v < 0.0)
            g.px[i] = 2.0 * v;
        else if (v > 1.0)
            g.px[i] = 2.0 * (v - 1.0);
        else
            g.px[i] = 0.0;
    }

    return g;
}

inline PenaltyResult penalty_method(
    const Image& raw,
    double lambda,
    double growth = 10.0,
    int max_outer = 7,
    int max_inner = 100,
    double tol = 1e-8
) {
    Image I = raw;

    double r = 1.0;
    double L = 1.0 + 8.0 * lambda;

    long calls = 0;
    long grad_calls = 0;

    std::vector<PenaltyRow> rows;

    for (int s = 0; s < max_outer; s++) {

        double Lr = L + 2.0 * r;
        double h = 1.0 / Lr;

        Image previous = I;
        int inner_iterations = 0;

        // Внутрішній метод Нестерова
        for (int k = 1; k <= max_inner; k++) {

            double beta = double(k - 1) / (k + 2);

            Image Y = I;

            for (size_t i = 0; i < I.px.size(); i++) {
                Y.px[i] = I.px[i]
                    + beta * (I.px[i] - previous.px[i]);
            }

            Image g = energy_grad(Y, raw, lambda);
            Image pg = penalty_grad(Y);

            for (size_t i = 0; i < g.px.size(); i++) {
                g.px[i] += r * pg.px[i];
            }

            grad_calls++;

            double gn = norm(g);

            if (gn <= 1e-6)
                break;

            Image next = Y;

            for (size_t i = 0; i < next.px.size(); i++) {
                next.px[i] -= h * g.px[i];
            }

            previous = I;
            I = next;

            inner_iterations++;
        }

        double E = energy(I, raw, lambda);
        double P = penalty_value(I);

        calls += 1;

        Violation v = violation(I);

        rows.push_back({
            s, r, Lr, h,
            E, P, r * P,
            v.max_value,
            v.count,
            inner_iterations,
            calls,
            grad_calls
        });

        if (r * P <= tol)
            break;

        r *= growth;
    }

    return {
        I,
        energy(I, raw, lambda),
        penalty_value(I),
        static_cast<int>(rows.size()),
        calls,
        grad_calls,
        rows
    };
}


// FISTA / accelerated projected gradient, based on the practicum's accelerated branch.
inline ProjectedResult fista_projected(const Image& raw, double lambda,
                                       int max_iter = 200, double tol = 1e-6) {
    double L = 1.0 + 8.0 * lambda;
    double h = 1.0 / L;
    Image I = project(raw), I_prev = I, Y = I;
    std::vector<ProjectedRow> rows;
    int calls = 0, grad_calls = 0;
    bool converged = false;
    double G = 0.0;
    for (int k=1; k<=max_iter; ++k) {
        Image grad = energy_grad(Y, raw, lambda);
        ++grad_calls;
        Image I_new = Y;
        for(size_t i=0;i<I_new.px.size();++i)
            I_new.px[i] = std::clamp(Y.px[i] - h*grad.px[i], 0.0, 1.0);
        double sum=0.0;
        for(size_t i=0;i<I_new.px.size();++i) {
            double d=(I_new.px[i]-Y.px[i])/h;
            sum+=d*d;
        }
        G = std::sqrt(sum);
        double E=energy(I,raw,lambda); ++calls;
        Violation v=violation(I);
        rows.push_back({k-1,E,G,active_constraints(I),v.max_value,calls,grad_calls});
        if (G <= tol) { converged=true; break; }
        double beta=double(k-1)/(k+2);
        for(size_t i=0;i<Y.px.size();++i)
            Y.px[i]=I_new.px[i]+beta*(I_new.px[i]-I_prev.px[i]);
        I_prev=I_new;
        I=I_new;
    }
    return {I, int(rows.size())-1, calls,grad_calls,
            energy(I,raw,lambda),G,converged,rows};
}
