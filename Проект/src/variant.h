
#pragma once

// Варіант №15
const int STUDENT_ID = 15;


const unsigned long long SEED =
    (STUDENT_ID * 7919ULL + 104729ULL) % 2147483647ULL;

// Параметри варіанта
const double NOISE_SIGMA = 0.1756;
const double X0 = 0.2492;
const double DELTA0 = 0.2866;
const double EPS = 0.0001;


const double LAMBDA0 = 0.6548;
const double SIMPLEX_STEP = 0.6420;


const int SIZE = 128;
