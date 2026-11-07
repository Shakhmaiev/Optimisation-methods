
#pragma once

#include <functional>

using namespace std;

struct Oracle {
    function<double(double)> f;
    int n_f = 0;

    Oracle(function<double(double)> func) {
        f = func;
    }

    // Звичайне обчислення з підрахунком
    double value(double x) {
        n_f++;
        return f(x);
    }

    // Обчислення без підрахунку
    double peek(double x) {
        return f(x);
    }

    void reset() {
        n_f = 0;
    }
};
