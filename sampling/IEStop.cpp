#pragma once
#include "IEStop.hpp"

double IEStop::mean() {
    double sum = 0;
    for (int i = 1; i <= t; i++) {
        sum += iearr[i];
    }
    sum /= t;
    return sum;
}

double IEStop::sd() {
    double ss = 0;
    double mean = IEStop::mean();
    for (int i = 1; i <= t; i++) {
        ss += (mean - iearr[i]) * (mean - iearr[i]);
    }
    ss /= t;
    return ss;
}

/*
execute one loop (from t-1 to t, do not change t)
return 0 if ready to stop
*/
int IEStop::loop() {
    double R = 1.0;
    if (t > floor(pow(BETA, k))) {
        k++;
        double alpha = floor(pow(BETA, k)) / floor(pow(BETA, k - 1));
        iex = -alpha * log((1 / pow(log(BETA, t), 1.1)) / 3)
    }
    double c = sd() * sqrt(2 * iex / t) + 3 * R * iex / t;
    lb = max(lb, mean() - c);
    ub = min(ub, mean() + c);
    if ((1 + eps) * lb >= (1 - eps) * ub) {
        return 0;
    }
    return 1;
}

void IEStop::init(double epsilon, double deltad) {
    memset(ieexi, 0, MAX_SAMPLE * sizeof(bool));
    memset(iearr, 0, MAX_SAMPLE * sizeof(int));
    t = 0;
    k = 0;
    lb = 0;
    ub = INFIN;
    eps = epsilon;
    delta = deltad;
}

/*
add X_rank to IEStop process
return 0 if stop
*/
int IEStop::add(int rank, int x) {
    iearr[rank] = x;
    ieexi[rank] = true;
    if (rank == t + 1) {
        t = rank;
        if (loop() == 0) {
            return 0;
        }
        while (ieexi[t + 1]) {
            t++;
            if (loop() == 0) {
                return 0;
            }
        }
    }
    return 1;
}

void IEStop::print_res() {
    double est = ((1 + eps) * lb + (1 - eps) * ub) / 2;
    printf("Estimation: %f.", est);
}