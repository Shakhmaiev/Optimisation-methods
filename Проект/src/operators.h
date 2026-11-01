
#pragma once

#include "image.h"
#include <cmath>
#include <cstddef>

struct Field {
    Image x;
    Image y;
};

// Дискретний градієнт
inline Field grad_xy(const Image& u) {

    Field g{
        Image(u.h, u.w),
        Image(u.h, u.w)
    };

    for (int i = 0; i < u.h; i++) {
        for (int j = 0; j < u.w; j++) {

            g.x.at(i, j) = (j < u.w - 1)
                ? u.at(i, j + 1) - u.at(i, j)
                : 0.0;

            g.y.at(i, j) = (i < u.h - 1)
                ? u.at(i + 1, j) - u.at(i, j)
                : 0.0;
        }
    }

    return g;
}

// Дивергенція
inline Image divergence(const Image& px,
                        const Image& py) {

    Image result(px.h, px.w);

    for (int i = 0; i < px.h; i++) {
        for (int j = 0; j < px.w; j++) {

            double value = 0.0;

            // Горизонтальний напрямок
            if (j == 0)
                value += px.at(i, 0);
            else if (j < px.w - 1)
                value += px.at(i, j) - px.at(i, j - 1);
            else
                value -= px.at(i, px.w - 2);

            // Вертикальний напрямок
            if (i == 0)
                value += py.at(0, j);
            else if (i < px.h - 1)
                value += py.at(i, j) - py.at(i - 1, j);
            else
                value -= py.at(px.h - 2, j);

            result.at(i, j) = value;
        }
    }

    return result;
}

// Лапласіан
inline Image laplacian(const Image& u) {

    Field g = grad_xy(u);

    return divergence(g.x, g.y);
}

// Енергія Тихонова
inline double energy(const Image& I,
                     const Image& noisy,
                     double lambda) {

    Field g = grad_xy(I);

    double fit = 0.0;
    double smooth = 0.0;

    for (size_t k = 0; k < I.px.size(); k++) {

        double d = I.px[k] - noisy.px[k];

        fit += d * d;

        smooth += g.x.px[k] * g.x.px[k]
                + g.y.px[k] * g.y.px[k];
    }

    return 0.5 * fit + 0.5 * lambda * smooth;
}

// Градієнт енергії Тихонова
inline Image energy_grad(const Image& I,
                         const Image& noisy,
                         double lambda) {

    Image lap = laplacian(I);
    Image result(I.h, I.w);

    for (size_t k = 0; k < I.px.size(); k++) {

        result.px[k] =
            (I.px[k] - noisy.px[k])
            - lambda * lap.px[k];
    }

    return result;
}

// Евклідова норма
inline double norm(const Image& I) {

    double sum = 0.0;

    for (double value : I.px)
        sum += value * value;

    return sqrt(sum);
}
