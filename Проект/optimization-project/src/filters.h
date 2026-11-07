
#pragma once

#include "image.h"
#include <vector>
#include <cmath>
#include <algorithm>

using namespace std;


Image gaussian_blur(const Image& img, double sigma) {

    if (sigma <= 0.000001)
        return img;

    // Радіус ядра
    int r = int(4.0 * sigma + 0.5);

    // Створення гаусового ядра
    vector<double> kernel(2 * r + 1);

    double sum = 0.0;

    for (int t = -r; t <= r; t++) {
        double value = exp(-(double(t) * t) /
                           (2.0 * sigma * sigma));

        kernel[t + r] = value;
        sum += value;
    }

    // Нормування ядра
    for (int t = 0; t < 2 * r + 1; t++) {
        kernel[t] /= sum;
    }

    Image temp(img.h, img.w);
    Image result(img.h, img.w);

    // горизонтальний
    for (int i = 0; i < img.h; i++) {
        for (int j = 0; j < img.w; j++) {

            double value = 0.0;

            for (int t = -r; t <= r; t++) {
                int jj = min(max(j + t, 0), img.w - 1);
                value += kernel[t + r] * img.at(i, jj);
            }

            temp.at(i, j) = value;
        }
    }

    // вертикальний
    for (int i = 0; i < img.h; i++) {
        for (int j = 0; j < img.w; j++) {

            double value = 0.0;

            for (int t = -r; t <= r; t++) {
                int ii = min(max(i + t, 0), img.h - 1);
                value += kernel[t + r] * temp.at(ii, j);
            }

            result.at(i, j) = value;
        }
    }

    return result;
}



// numpy.pad(..., mode="reflect") не дублює крайній піксель.
// -1 -> 1, size -> size-2. Застосовується лише в bilateral_filter.
inline int mirror(int index, int size) {
    if (size <= 1) return 0;
    while (index < 0 || index >= size) {
        if (index < 0) index = -index;
        if (index >= size) index = 2 * size - index - 2;
    }
    return index;
}

// Білатеральний фільтр
inline Image bilateral_filter(const Image& img,
                              double sigma_spatial,
                              double sigma_range) {

    sigma_spatial = max(sigma_spatial, 0.001);
    sigma_range = max(sigma_range, 0.001);

    // Радіус вікна
    int r = clamp(
        int(lround(2.0 * sigma_spatial)), 1, 5
    );

    Image result(img.h, img.w);

    double two_ss = 2.0 * sigma_spatial * sigma_spatial;
    double two_sr = 2.0 * sigma_range * sigma_range;

    for (int i = 0; i < img.h; i++) {
        for (int j = 0; j < img.w; j++) {

            double center = img.at(i, j);

            double numerator = 0.0;
            double denominator = 0.0;

            // Перебираємо сусідні пікселі
            for (int dy = -r; dy <= r; dy++) {
                for (int dx = -r; dx <= r; dx++) {

                    int y = mirror(i + dy, img.h);
                    int x = mirror(j + dx, img.w);

                    double neighbor = img.at(y, x);

                    // Просторова вага
                    double ws = exp(
                        -(dx * dx + dy * dy) / two_ss
                    );

                    // Вага за яскравістю
                    double diff = neighbor - center;

                    double wr = exp(
                        -(diff * diff) / two_sr
                    );

                    double weight = ws * wr;

                    numerator += weight * neighbor;
                    denominator += weight;
                }
            }

            result.at(i, j) = numerator / denominator;
        }
    }

    return result;
}
