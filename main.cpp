#include <iostream>
#include <chrono>
#include <vector>
#include <fstream>
#include "instance.h"
#include "l1.h"
#include "cutting_plane.h"
//#include "omp.h"
//#include "DCA.h"


int main() {
    using namespace chrono;
    Instance inst(100, 2000);

    system_clock::time_point t1 = high_resolution_clock::now();


    system_clock::time_point t2 = high_resolution_clock::now();
    duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
//    cout << "生成随机实例用时：
    vector<int> success_frequency_l1;
    vector<int> success_frequency_weighted_l1;
    vector<int> success_frequency_cp;

    for (int sparse_level = 5; sparse_level < 40; sparse_level += 5) {
        cout << sparse_level << endl;
        int n_success_l1 = 0;
        int n_success_weighted_l1 = 0;
        int n_success_cp = 0;
        for (int i = 0; i < 100; ++i) {
            generate_instance(inst, sparse_level);
            VectorXd x = inst.x;

            inst.clear_solution();
            l1_method(inst);
            if ((inst.x - x).norm() / x.norm() <= 1e-2) {
                n_success_l1 += 1;
            }

            inst.clear_solution();
            reweighted_l1(inst);
            if ((inst.x - x).norm() / x.norm() <= 1e-2) {
                n_success_weighted_l1 += 1;
            }

            inst.clear_solution();
            cutting_plane_method(inst);
            if ((inst.x - x).norm() / x.norm() <= 1e-2) {
                n_success_cp += 1;
            }
        }
        success_frequency_l1.emplace_back(n_success_l1);
        success_frequency_weighted_l1.emplace_back(n_success_weighted_l1);
        success_frequency_cp.emplace_back(n_success_cp);
    }

    for (auto &x:success_frequency_l1) {
        cout << x << ",";
    }
    cout << endl;
    for (auto &x:success_frequency_weighted_l1) {
        cout << x << ",";
    }
    cout << endl;
    for (auto &x:success_frequency_cp) {
        cout << x << ",";
    }
    cout << endl;
//    t1 = high_resolution_clock::now();
//    l1_method(inst);
//    t2 = high_resolution_clock::now();
//    time_span = duration_cast<duration<double>>(t2 - t1);
//    cout << time_span.count() << "s" << endl;
//    cout << inst.count_nnz() << endl;
//
//
//    inst.clear_solution();
//    t1 = high_resolution_clock::now();
//    reweighted_l1(inst);
//    t2 = high_resolution_clock::now();
//    time_span = duration_cast<duration<double>>(t2 - t1);
//    cout << time_span.count() << "s" << endl;
//    cout << inst.count_nnz() << endl;

    return 0;
}