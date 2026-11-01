#pragma once

#include <vector>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <stdexcept>
#include <string>

#include "rng.h"

using namespace std;


// Структура зображення
struct Image {
    int h = 0, w = 0;
    vector<double> px;

    Image() = default;

    Image(int height, int width, double value = 0.0) {
        h = height;
        w = width;
        px.resize(h * w, value);
    }

    double& at(int i, int j) {
        return px[i * w + j];
    }

    const double& at(int i, int j) const {
        return px[i * w + j];
    }
};


// Еталонне зображення
Image clean_image(int size = 128) {
    Image img(size, size);

    const double PI = 3.14159265358979323846;

    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {

            double u = double(j) / (size - 1);
            double v = double(i) / (size - 1);

            // 1. Фон
            double val = 0.25 + 0.35 * (0.6 * u + 0.4 * v);
            val += 0.08 * cos(2.0 * PI *
                sqrt((u - 0.5) * (u - 0.5) +
                     (v - 0.5) * (v - 0.5)));

            // 2. Світле коло
            if ((u - 0.32) * (u - 0.32) +
                (v - 0.34) * (v - 0.34) <= 0.16 * 0.16) {
                val = 0.92;
            }

            // 3. Темний квадрат
            if (abs(u - 0.70) <= 0.15 &&
                abs(v - 0.30) <= 0.15) {
                val = 0.10;
            }

            // 4. Трикутник
            if (v > 0.60 && v < 0.90 &&
                abs(u - 0.30) < 0.9 * (v - 0.60)) {
                val = 0.68;
            }

            // 5. Шахівниця
            if (u > 0.58 && u < 0.92 &&
                v > 0.60 && v < 0.90) {

                int cell = max(3, size / 12);

                if (((j / cell) + (i / cell)) % 2 == 0)
                    val = 0.85;
                else
                    val = 0.20;
            }

            img.at(i, j) = clamp(val, 0.0, 1.0);
        }
    }

    // 6. Чотири білі смуги вгорі
    for (int idx = 0; idx < 4; idx++) {
        int thickness = max(2, (4 - idx) * max(1, size / 64));
        int row = int(0.04 * size) + idx * max(3, size / 22);

        if (row + thickness < size) {
            for (int i = row; i < row + thickness; i++) {
                for (int j = int(0.05 * size); j < int(0.95 * size); j++) {
                    img.at(i, j) = 0.98;
                }
            }
        }
    }

    return img;
}


// Зашумлене зображення
// raw  - без обрізання
// noisy - після обрізання до [0,1]

struct NoisyData {
    Image raw;
    Image noisy;
};

NoisyData add_noise(const Image& clean, double sigma, uint64_t seed) {
    Rng rng(seed);

    Image raw(clean.h, clean.w);
    Image noisy(clean.h, clean.w);

    for (int i = 0; i < clean.h; i++) {
        for (int j = 0; j < clean.w; j++) {
            double value = clean.at(i, j) + sigma * rng.normal();

            raw.at(i, j) = value;
            noisy.at(i, j) = clamp(value, 0.0, 1.0);
        }
    }

    return {raw, noisy};
}


// Збереження у PGM
void save_pgm(const Image& img, const string& filename) {
    ofstream file(filename, ios::binary);

    if (!file) {
        throw runtime_error("Cannot open file: " + filename);
    }

    file << "P5\n" << img.w << " " << img.h << "\n255\n";

    for (double x : img.px) {
        unsigned char b = (unsigned char)lround(clamp(x, 0.0, 1.0) * 255.0);
        file.write((const char*)&b, 1);
    }
}