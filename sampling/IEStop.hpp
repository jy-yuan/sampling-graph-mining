#ifndef IESTOP_HPP
#define IESTOP_HPP
#pragma once
#include <cmath>
#include <cstdio>
#include <cstring>

#define MAX_SAMPLE 10000
#define BETA 10
#define INFIN 1 << 20

/*
Interval Estimation Stopping
based on central limit theorem
*/

class IEStop {
   public:
    static IEStop& get_instance() {
        static IEStop ieins;
        return ieins;
    }
    void init(double epsilon, double deltad);
    int add(int rank, int x);
    void print_res();

   private:
    bool ieexi[MAX_SAMPLE];
    int iearr[MAX_SAMPLE];
    int t = 0;
    int k = 0;
    double lb = 0;
    double ub = INFIN;
    double eps = 0.05;
    double delta = 0.05;
    double iex = 1;
    IEStop();
    ~IEStop(){};
    double mean();
    double sd();
    int loop();
};

#endif