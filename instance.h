//
// Created by agamemnon on 17-9-1.
//

#ifndef SIGNAL_RECOVERY_INSTANCE_H
#define SIGNAL_RECOVERY_INSTANCE_H

#include <iostream>
#include <Eigen/Dense>
#include <random>

using namespace std;
using namespace Eigen;

struct Instance {
    MatrixXd A;
    VectorXd b;
    VectorXd x;

    Instance(const int m, const int n) {
        A = MatrixXd::Zero(m, n);
        b = VectorXd::Zero(m);
        x = VectorXd::Zero(n);
    }

    void display() {
        cout << "A" << endl;
        cout << A << endl;
        cout << "b" << endl;
        cout << b.transpose() << endl;
    }

    int count_nnz() {
        int nnz = 0;
        for (int j = 0; j < x.rows(); ++j) {
            if (abs(x[j]) > 1e-4) {
                nnz++;
            }
        }
        return nnz;
    }

    bool check_feasiblity() {
        VectorXd y = A * x;
        bool flag = (y - b).norm() <= 1e-4;
        return flag;
    }

    void clear_solution() {
        x = VectorXd::Zero(x.rows());
    }

    void initialize_solution(VectorXd &x0) {
        x = x0;
    }
};


void generate_instance(Instance &inst, int sparse_level, double miu = 0, double sigma = 1) {
    random_device rd;
    default_random_engine dre(rd());
    int m = inst.A.rows();
    int n = inst.A.cols();
    normal_distribution<double> nd(miu, sigma / n);
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            inst.A(i, j) = nd(dre);
        }
    }

    VectorXd x = VectorXd::Zero(n);
    int nnz = 0;
    uniform_int_distribution<int> uid(0, n - 1);
    normal_distribution<double> nd1(miu, sigma);
    while (nnz < sparse_level) {
        int j = uid(dre);
        if (x(j) == 0) {
            x(j) = nd1(dre);
            nnz += 1;
        }
    }

    inst.b = inst.A * x;
    inst.x = x;
}

#endif //SIGNAL_RECOVERY_INSTANCE_H
