#pragma once
#include <bitset>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <functional>
#include <map>
#include <tuple>
#include <vector>

#include "logger.hpp"

struct tableRow {
    /// @todo write the constructor for init
    // one row of table for nm output (see actual task)
    tableRow(double xi,
             double vi,
             double yi,
             double v2i,
             double y2i,
             double viv2i,
             double yiy2i,
             double LE,
             double LEOLE,
             double hi,
             uint C1,
             uint C2,
             double ui,
             double uvi)
        : xi(xi), vi(vi), yi(yi), v2i(v2i), y2i(y2i), viv2i(viv2i), yiy2i(yiy2i), LE(LE), LEOLE(LEOLE), hi(hi), C1(C1), C2(C2), ui(ui), uvi(uvi) {}
    double xi;
    double vi;
    double yi;     // only nedded if there is a system of two variables
    double v2i;    // only needed if LEC is on
    double y2i;    // only nedded if there is a system of two variables and lec is on
    double viv2i;  // only needed if LEC is on and
    double yiy2i;  // only nedded if there is a system of two variables and lec is on
    double LE;     // only visible in table if LEC is on
    double LEOLE;  // same
    double hi;
    uint C1;     // this values actually is constant when lec is off
    uint C2;     //
    double ui;   // This cols visible only if there is analytic solution in task (so set it invisible else)
    double uvi;  //
    // 14 cols

    std::tuple<double, double, double, double, double, double, double, double, double, double, uint, uint, double, double> get_tuple() {
        return std::tuple<double, double, double, double, double, double, double, double, double, double, uint, uint, double, double>(xi,
                                                                                                                                      vi,
                                                                                                                                      yi,
                                                                                                                                      v2i,
                                                                                                                                      y2i,
                                                                                                                                      viv2i,
                                                                                                                                      yiy2i,
                                                                                                                                      LE,
                                                                                                                                      LEOLE,
                                                                                                                                      hi,
                                                                                                                                      C1,
                                                                                                                                      C2,
                                                                                                                                      ui,
                                                                                                                                      uvi);
    }

    friend std::ostream& operator<<(std::ostream& os, const tableRow& row) {
        os << "xi = " << row.xi << ", vi = " << row.vi << ", yi = " << row.yi << ", y2i = " << row.y2i << ", v2i = " << row.v2i << ", viv2i = " << row.viv2i
           << ", LE = " << row.LE << ", hi = " << row.hi << ", C1 = " << row.C1 << ", C2 = " << row.C2 << ", ui = " << row.ui << ", uvi = " << row.uvi;
        return os;
    }

    friend std::ofstream& operator<<(std::ofstream& os, const tableRow& row) {
        os << row.xi << '\t' << row.vi << '\t' << row.yi << '\t' << row.y2i << '\t' << row.v2i << '\t' << row.viv2i << '\t' << row.LE << '\t' << row.hi << '\t'
           << row.C1 << '\t' << row.C2 << '\t' << row.ui << '\t' << row.uvi;
        return os;
    }
};

typedef std::vector<tableRow> resultTable;  // table for output

struct config {
    config() = default;
    config(double x_min,
           double x_max,
           double x_0,
           double u_0,
           double du_0,
           double step,
           uint N_max,
           bool LEC,
           double eps,
           double A,
           double B,
           double C,
           double brd_prec,
           std::bitset<14> isColVis)
        : x_min(x_min), x_max(x_max), x_0(x_0), u_0(u_0), du_0(du_0), step(step), N_max(N_max), LEC(LEC), eps(eps), A(A), B(B), C(C), brdr_prec(brd_prec), isColVisible(isColVis) {}
    config(std::tuple<double, double, double, double, double, double, uint, bool, double, double, double, double, double> tpl)
        : config(std::get<0>(tpl),
                 std::get<1>(tpl),
                 std::get<2>(tpl),
                 std::get<3>(tpl),
                 std::get<4>(tpl),
                 std::get<5>(tpl),
                 std::get<6>(tpl),
                 std::get<7>(tpl),
                 std::get<8>(tpl),
                 std::get<9>(tpl),
                 std::get<10>(tpl),
                 std::get<11>(tpl),
                 std::get<12>(tpl),
                 0b11111111111111) {}
    // Left and right limits for x variable
    double x_min = 0.l;
    double x_max = 0.l;

    // initial point
    double x_0 = 0.l;
    double u_0 = 0.l;
    double du_0 = 0.l;

    double step = 0.l;  // step
    uint N_max = 0;   // Maximum num for iterations

    bool LEC = 1;      // Is there control for local error
    double eps = 0.f;  // Epsilon for local error control

    double A = 0.l;
    double B = 0.l;
    double C = 0.l;

    double brdr_prec = 0.l; // right boarder precision

    std::bitset<14> isColVisible = { 0b11111111111111 };  // bit set which keeps visibility of each column in table (0 is not visible and 1 is visible)
    // keep in mind that order of bits is inverted!!!
    bool is_col_visible(const size_t& col) const {
        return isColVisible[col];
    }

    void set_col_visible(const size_t& col) {
        isColVisible[col] = true;
    }

    void set_col_invisible(const size_t& col) {
        isColVisible[col] = false;
    }

    constexpr void set_col_visibility(std::bitset<14> val) {
        isColVisible = val;
    }

    friend std::ostream& operator<<(std::ostream& os, const config& cfg) {
        os << "x_min=" << cfg.x_min << " x_max=" << cfg.x_max << " x_0=" << cfg.x_0 << " u_0=" << cfg.u_0 << " step=" << cfg.step << " N_max=" << cfg.N_max
           << " LEC=" << cfg.LEC << " eps=" << cfg.eps;
        return os;
    }
};

namespace utils {
// /// @brief make config from python input
// config make_config(const double& x_min,
//                    const double& x_max,
//                    const double& x_0,
//                    const double& u_0,
//                    const double& du_0,
//                    const double& step,
//                    const uint& N_max,
//                    const bool& LEC,
//                    const double& eps,
//                    const double& A,
//                    const double& B,
//                    const double& C,
//                    std::bitset<14>);

/// numerical method
resultTable RK4(std::function<double(double, double)> rhs, const config& cfg);
resultTable RK3(std::function<double(double, double)>&& rhs, const config& cfg);

resultTable RK3_SOE(std::function<double(double, double, double)>&& rhs1, std::function<double(double, double, double)>&& rhs2, const config& cfg);

/// @brief  RK4 for System Of two Equations
resultTable RK4_SOE(std::function<double(double, double, double)> rhs1, std::function<double(double, double, double)> rhs2, const config& cfg);

/// step for method
inline double StepRK4(std::function<double(double, double)> rhs, const double& x, const double& u, const double& step);

/// step RK4 for System Of two Equations
inline std::vector<double> StepRK4_SOE(std::function<double(double, double, double)> rhs1,
                                       std::function<double(double, double, double)> rhs2,
                                       const double& x,
                                       const double& u,
                                       const double& y,
                                       const double& step);

resultTable RK4_LS(std::function<double(double, double, double)> rhs, const config& cfg);

}  // namespace utils

template <class Func, class Tuple, std::size_t... Is>
void apply_elemwise(Func&& f, const Tuple& tp, std::index_sequence<Is...>) {
    (f(std::get<Is>(tp)), ...);
}

static std::function<resultTable(std::function<double(double, double)>, config)> task_rk4 = utils::RK4;

static std::function<resultTable(std::function<double(double, double, double)>, config)> task_rk4_lseq = utils::RK4_LS;

std::function<double(const double&)> make_test_true_sol(const double& x_0, const double& u_0);

void calculate_global_error(resultTable& tbl, std::function<double(const double&)>&& f);

std::tuple<double, double> find_max_LE(const resultTable& tbl);

std::tuple<double, double> find_max_h(const resultTable& tbl);

std::tuple<double, double> find_min_h(const resultTable& tbl);

std::tuple<double, double> find_max_uvi(const resultTable& tbl);

void dumpTableInFile(const resultTable&, std::ofstream&);