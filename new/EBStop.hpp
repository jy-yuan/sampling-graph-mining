#pragma once
#include <cstdio>
#include <cstring>

/*
Empirical Bernstein Stopping
@inproceedings{mnih2008empirical,
  title={Empirical bernstein stopping},
  author={Mnih, Volodymyr and Szepesv{\'a}ri, Csaba and Audibert, Jean-Yves},
  booktitle={Proceedings of the 25th international conference on Machine learning},
  pages={672--679},
  year={2008}
}
*/

class EBStop {
  public:
    static EBStop& get_instance() {
        static EBStop ebins;
        return ebins;
    }
    void init();
    int add(int rank, int x);
    void print_res();
  private:
    bool ebexi[100];
    int ebarr[100];
    EBStop();
    ~EBStop() {};
};