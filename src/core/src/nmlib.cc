#include "nmlib.hpp"
#include <bitset>

// config utils::make_config(const double& x_min,
//                           const double& x_max,
//                           const double& x_0,
//                           const double& u_0,
//                           const double& du_0,
//                           const double& step,
//                           const uint& N_max,
//                           const bool& LEC,
//                           const double& eps,
//                           const double& A,
//                           const double& B,
//                           const double& C,
//                           std::bitset<14> cols) {
//     LOG_INFO_CLI("Making some config");
// #if defined(DEBUG)
//     LOG_DEBUG_CLI("Making config from python data", "Checking arguments");
// #endif
//     NM_ASSERT((x_max >= x_min), "invalid argument x_min > x_max");
//     /// @todo write an asserts for other args

//     config cfg{ x_min, x_max, x_0, u_0, du_0, step, N_max, LEC, eps, A, B, C, cols };
//     LOG_INFO_CLI("Configuration done");
//     return cfg;
// }

std::function<double(const double&)> make_test_true_sol(const double& x_0, const double& u_0) {
    return [&](const double& x) {
        return u_0 * std::exp(2 * (x - x_0));
    };
}

void calculate_global_error(resultTable& tbl, std::function<double(const double&)>&& f) {
    for (auto& row : tbl) {
        row.ui = f(row.xi);
        row.uvi = std::abs(row.vi - row.ui);
    }
}

std::tuple<double, double> find_max_LE(const resultTable& tbl) {
    double max = 0.l;
    double x = 0.l;

    for (auto& row : tbl) {
        double tmp = std::abs(row.LE);
        if (tmp > max){
            max = tmp;
            x = row.xi;
        }
    }

    return {max, x};
}

std::tuple<double, double> find_max_h(const resultTable& tbl) {
    double max = 0.f;
    double x = 0.l;

    for (auto& row : tbl) {
        if (row.hi > max) {
            max = row.hi;
            x = row.xi;
        }
    }

    return {max, x};
}

std::tuple<double, double> find_min_h(const resultTable& tbl) {
    double min = tbl.at(0).hi;
    double x = 0.l;

    for (auto& row : tbl) {
        if (row.hi < min) {
            min = row.hi;
            x = row.xi;
        }
    }

    return {min, x};
}

std::tuple<double, double> find_max_uvi(const resultTable& tbl) {
    double max = 0.f;
    double x = 0.l;

    for (auto& row : tbl) {
        if (row.uvi > max) {
            max = row.uvi;
            x = row.xi;
        }
    }

    return {max, x};
}

void dumpTableInFile(const resultTable& tbl, std::ofstream& fs) {
    for (const auto& row : tbl) {
        fs << row << '\n';
    }
}