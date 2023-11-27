#include "nmlib.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <tuple>
#include <utility>

inline double utils::StepRK4(std::function<double(double, double)> rhs, const double& x, const double& u, const double& step) {
    LOG_INFO_CLI("Start RK4 step with following config", x, u, step);

    double new_u, k1, k2, k3, k4;
    k1 = rhs(x, u);
    k2 = rhs(x + step / 2, u + step / 2 * k1);
    k3 = rhs(x + step / 2, u + step / 2 * k2);
    k4 = rhs(x + step, u + step * k3);
    new_u = u + step / 6 * (k1 + (2 * k2) + (2 * k3) + k4);
    return (new_u);
}

resultTable utils::RK4(std::function<double(double, double)> rhs, const config& cfg) {
    LOG_INFO_CLI("Start RK4 with following config", cfg);

    resultTable table;
    double xi, x_min, x_max, ui, stepi, N_max, eps;
    uint C1, C2;
    C1 = 0;
    C2 = 0;
    xi = cfg.x_0;
    x_min = cfg.x_min;
    x_max = cfg.x_max;
    ui = cfg.u_0;
    stepi = cfg.step;
    N_max = cfg.N_max;
    eps = cfg.eps;

    if (cfg.LEC) {
        int i = 0;
        double u1, u2;
        double viv2i, LocalError, Old_LocalError;
        while (xi >= x_min && xi + stepi < x_max && i <= N_max) {
            u1 = StepRK4(rhs, xi, ui, stepi);
            u2 = StepRK4(rhs, xi, ui, stepi / 2);
            u2 = StepRK4(rhs, xi + stepi / 2, u2, stepi / 2);
            Old_LocalError = LocalError;
            LocalError = (std::abs(u1 - u2)) / (std::pow(2, 4) - 1);
            if (LocalError < eps / std::pow(2, 5)) {
                ui = u1;
                xi = xi + stepi;
                viv2i = std::abs(ui - u2);
                tableRow row(xi, ui, 0.f, u2, viv2i, 0.f, 0.f, LocalError, LocalError / Old_LocalError, stepi, C1, C2, 0.f, 0.f);
                table.push_back(row);
                stepi = stepi * 2;
                C2++;
                i++;
            } else if (LocalError >= eps / std::pow(2, 5) && LocalError <= eps) {
                ui = u1;
                xi = xi + stepi;
                viv2i = std::abs(ui - u2);
                tableRow row(xi, ui, 0.f, u2, 0.f, viv2i, 0.f, LocalError, LocalError / Old_LocalError, stepi, C1, C2, 0.f, 0.f);
                table.push_back(row);
                i++;

            } else if (LocalError > eps) {
                stepi = stepi / 2;
                C1++;
                i++;
            } else {
                LOG_ERROR_CLI(cfg);
            }
        }
        LOG_INFO_CLI("RK4 close", cfg);
        return (table);

    } else if (not(cfg.LEC)) {
        int i = 0;
        while (xi >= x_min && xi + stepi < x_max && i <= N_max) {
            ui = StepRK4(rhs, xi, ui, stepi);
            xi = xi + stepi;
            i++;
            tableRow row(xi, ui, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, stepi, 0.f, 0.f, 0.f, 0.f);
            table.push_back(row);
        }
        LOG_INFO_CLI("RK4 close", cfg);
        return (table);
    } else {
        LOG_INFO_CLI("Error in RK4", cfg);
    }
}

inline std::vector<double> utils::StepRK4_SOE(std::function<double(double, double, double)> rhs1,
                                              std::function<double(double, double, double)> rhs2,
                                              const double& x,
                                              const double& u,
                                              const double& y,
                                              const double& step) {
    LOG_DEBUG_CLI("Start StepRK4_SOE step with following config", x, u, step);

    double new_u, new_y;
    std::vector<double> res;
    new_u = u + step * rhs1(x + step / 2, u + (step / 2) * (rhs1(x, u, y)), y + (step / 2) * (rhs2(x, u, y)));
    new_y = y + step * rhs2(x + step / 2, u + (step / 2) * (rhs1(x, u, y)), y + (step / 2) * (rhs2(x, u, y)));
    res.push_back(new_u);
    res.push_back(new_y);
    return (res);
}

resultTable utils::RK4_SOE(std::function<double(double, double, double)> rhs1, std::function<double(double, double, double)> rhs2, const config& cfg) {
    LOG_INFO_CLI("Start RK4_SOE with following config", cfg);
    resultTable table;
    double xi, x_min, x_max, ui, yi, stepi, N_max, eps;
    uint C1, C2;
    std::vector<double> tmp1, tmp2;
    C1 = 0;
    C2 = 0;
    xi = cfg.x_0;
    x_min = cfg.x_min;
    x_max = cfg.x_max;
    double a = cfg.A;
    double b = cfg.B;
    double c = cfg.C;
    ui = cfg.u_0;
    yi = cfg.du_0;
    stepi = cfg.step;
    N_max = cfg.N_max;
    eps = cfg.eps;
    tableRow row(xi, ui, yi, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, C1, C2, 0.f, 0.f);
    table.push_back(row);
    if (cfg.LEC) {
        double LocalError, Old_LocalError, viv2i, yiy2i;
        int i = 0;
        while (xi >= x_min && xi + stepi < x_max && i <= N_max) {
            //LOG_DEBUG_CLI("Start RK4_SOE with localstecpcontrol", cfg);
            tmp1 = StepRK4_SOE(rhs1, rhs2, xi, ui, yi, stepi);
            tmp2 = StepRK4_SOE(rhs1, rhs2, xi, ui, yi, stepi / 2);
            tmp2 = StepRK4_SOE(rhs1, rhs2, xi, tmp2.at(0), tmp2.at(1), stepi / 2);
            Old_LocalError = LocalError;
            LocalError = std::max(std::abs(tmp1.at(0) - tmp2.at(0)), std::abs(tmp1.at(1) - tmp2.at(1))) / 3;
            //LOG_DEBUG_CLI("tmp generated", cfg);
            if (LocalError < eps / std::pow(2, 5)) {
                //LOG_DEBUG_CLI("try LocalError < lec", cfg);
                yi = tmp1.at(1);
                ui = tmp1.at(0);
                xi = xi + stepi;
                viv2i = std::abs(ui - tmp2.at(0));
                yiy2i = std::abs(yi - tmp2.at(1));
                //LOG_DEBUG_CLI("viv2i generated", cfg);
                tableRow row(xi,
                             ui,
                             yi,
                             tmp1.at(0),
                             tmp2.at(1),
                             viv2i,
                             yiy2i,
                             LocalError,
                             LocalError / Old_LocalError,
                             stepi,
                             C1,
                             C2,
                             std::abs((10 / std::exp(xi / 100) - 3 / std::exp(1000 * xi)) - ui),
                             std::abs((3 / std::exp(1000 * xi) + 10 / std::exp(xi / 100)) - yi));
                //LOG_DEBUG_CLI("try push_back", cfg);
                table.push_back(row);
                //LOG_DEBUG_CLI("push_back done", cfg);
                stepi = stepi * 2;
                C2++;
                i++;
                //LOG_DEBUG_CLI("LocalError < lec done", cfg);
            } else if (LocalError >= eps / std::pow(2, 5) && LocalError <= eps) {
                //LOG_DEBUG_CLI("LocalError >= lec", cfg);
                yi = tmp1.at(1);
                ui = tmp1.at(0);
                xi = xi + stepi;
                viv2i = std::abs(ui - tmp2.at(0));
                yiy2i = std::abs(yi - tmp2.at(1));
                tableRow row(xi,
                             ui,
                             yi,
                             tmp2.at(0),
                             tmp2.at(1),
                             viv2i,
                             yiy2i,
                             LocalError * 4,
                             LocalError / Old_LocalError,
                             stepi,
                             C1,
                             C2,
                             std::abs((10 / std::exp(xi / 100) - 3 / std::exp(1000 * xi)) - ui),
                             std::abs((3 / std::exp(1000 * xi) + 10 / std::exp(xi / 100)) - yi));
                table.push_back(row);
                i++;
            } else if (LocalError > eps) {
                LOG_DEBUG_CLI("LocalError > lec*", cfg);
                stepi = stepi / 2;
                C1++;
                i++;
            } else {
                LOG_ERROR_CLI(cfg);
            }
        }
        return (table);
    } else if (not(cfg.LEC)) {
        //LOG_DEBUG_CLI("Start RK4_SOE without localstecpcontrol", cfg);
        int i = 0;
        while (xi >= x_min && xi + stepi < x_max && i <= N_max) {
            tmp1 = StepRK4_SOE(rhs1, rhs2, xi, ui, yi, stepi);
            ui = tmp1.at(0);
            yi = tmp1.at(1);
            xi = xi + stepi;
            i++;
            tableRow row(xi,
                         ui,
                         yi,
                         0.f,
                         0.f,
                         0.f,
                         0.f,
                         0.f,
                         0.f,
                         stepi,
                         0.f,
                         0.f,
                         std::abs((10 / std::exp(xi / 100) - 3 / std::exp(1000 * xi)) - ui),
                         std::abs((3 / std::exp(1000 * xi) + 10 / std::exp(xi / 100)) - yi));
            table.push_back(row);
        }
        return (table);
    } else {
        LOG_INFO_CLI("Error in RK4", cfg);
    }
}

resultTable utils::RK4_LS(std::function<double(double, double, double)> rhs, const config& cfg) {
    return utils::RK4_SOE(
        rhs,
        [&](double x, double u, double du) {
            return (499.995 * du - 500.005 * u);
        },
        cfg);
}

inline double RK3_STEP(std::function<double(double, double)>&& rhs, const double& step, const double& x, const double& v) {
    /* 
    x_(n+1) = x_n + h_n,
    v_(n+1) = v_n + h_n/6 * (k1 + 4k2 + k3),
    k1 = f(x_n, v_n),
    k2 = f(x_n + h_n/2, v_n + h_n/2 * k1),
    k3 = f(x_n + h_n, v_n + h_n * (-k1 + 2k2))
    */
    double k1 = rhs(x, v);
    double k2 = rhs(x + step / 2, v + step / 2 * k1);
    double k3 = rhs(x + step, v + step * (-k1 + 2 * k2));

    return v + step / 6 * (k1 + 4 * k2 + k3);
}

resultTable utils::RK3(std::function<double(double, double)>&& rhs, const config& cfg) {
    resultTable table;
    double xi, x_min, x_max, vi, stepi, N_max, eps;
    uint C1 = 0, C2 = 0;
    xi = cfg.x_0;
    x_min = cfg.x_min;
    x_max = cfg.x_max;
    vi = cfg.u_0;
    stepi = cfg.step;
    N_max = cfg.N_max;
    eps = cfg.eps;

    if (!cfg.LEC) {
        size_t i = 0;
        while (xi >= x_min && xi + stepi < x_max && i <= N_max) {
            vi = RK3_STEP(std::move(rhs), stepi, xi, vi);
            xi = xi + stepi;
            table.push_back({ xi, vi, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, stepi, 0, 0, 0.f, 0.f });
        }
        return table;
    } else {
        int i = 0;
        double v1, v2;
        double viv2i, LocalError = 0.f, Old_LocalError;
        while (xi >= x_min && xi + stepi < x_max && i <= N_max) {
            v1 = RK3_STEP(std::move(rhs), stepi, xi, vi);
            v2 = RK3_STEP(std::move(rhs), stepi / 2, xi, vi);
            v2 = RK3_STEP(std::move(rhs), stepi / 2, xi + stepi / 2, v2);
            Old_LocalError = LocalError;
            LocalError = (std::abs(v1 - v2)) / 7;
            if (LocalError < eps / 16) {
                vi = v1;
                xi = xi + stepi;
                viv2i = std::abs(vi - v2);
                table.push_back({ xi, vi, 0.f, v2, viv2i, 0.f, 0.f, LocalError, LocalError / Old_LocalError, stepi, C1, C2, 0.f, 0.f });
                stepi = stepi * 2;
                C2++;
                i++;
            } else if (LocalError >= eps / 16 && LocalError <= eps) {
                vi = v1;
                xi = xi + stepi;
                viv2i = std::abs(vi - v2);
                table.push_back({ xi, vi, 0.f, v2, 0.f, viv2i, 0.f, LocalError, LocalError / Old_LocalError, stepi, C1, C2, 0.f, 0.f });
                i++;

            } else if (LocalError > eps) {
                stepi = stepi / 2;
                C1++;
                i++;
            }
            return table;
        }
    }
    return table;
}

inline std::tuple<double, double> RK3_SOE_STEP(std::function<double(double, double, double)>&& rhs1,
                                               std::function<double(double, double, double)>&& rhs2,
                                               const double& step,
                                               const double& x,
                                               const double& vi1,
                                               const double& vi2) {
    /* 
    x_(n+1) = x_n + h_n,
    v_(n+1) = v_n + h_n/6 * (k1 + 4k2 + k3),
    k1 = f(x_n, v_n),
    k2 = f(x_n + h_n/2, v_n + h_n/2 * k1),
    k3 = f(x_n + h_n, v_n + h_n * (-k1 + 2k2))
    */
    double vi1_new = vi1 + step / 6 *
                               (rhs1(x, vi1, vi2) + 4 * (rhs1(x + step / 2, vi1 + step / 2 * rhs1(x, vi1, vi2), vi2)) +
                                rhs1(x + step, vi1 + step * (-rhs1(x, vi1, vi2) + 2 * (rhs1(x + step / 2, vi1 + step / 2 * rhs1(x, vi1, step), vi2))), vi2));
    double vi2_new = vi2 + step / 6 *
                               (rhs2(x, vi1, vi2) + 4 * (rhs2(x + step / 2, vi1, vi2 + step / 2 * rhs2(x, vi1, vi2))) +
                                rhs2(x + step, vi1, vi2 + +step * (-rhs2(x, vi1, vi2) + 2 * (rhs2(x + step / 2, vi1, vi2 + step / 2 * rhs2(x, vi1, step))))));

    return { vi1_new, vi2_new };
}

resultTable utils::RK3_SOE(std::function<double(double, double, double)>&& rhs1, std::function<double(double, double, double)>&& rhs2, const config& cfg) {
    resultTable table;
    double xi, x_min, x_max, vi, yi, stepi, N_max, eps;
    uint C1 = 0, C2 = 0;
    xi = cfg.x_0;
    x_min = cfg.x_min;
    x_max = cfg.x_max;
    vi = cfg.u_0;
    yi = cfg.du_0;
    stepi = cfg.step;
    N_max = cfg.N_max;
    eps = cfg.eps;

    if (!cfg.LEC) {
        size_t i = 0;
        while (xi >= x_min && xi + stepi < x_max && i <= N_max) {
            auto [vi_new, yi_new] = RK3_SOE_STEP(std::move(rhs1), std::move(rhs2), stepi, xi, vi, yi);
            xi = xi + stepi;
            vi = vi_new;
            yi = yi_new;

            table.push_back({ xi, vi, 0.f, yi, 0.f, 0.f, 0.f, 0.f, 0.f, stepi, 0, 0, 0.f, 0.f });
        }
        return table;
    } else {
        double LocalError, Old_LocalError, viv2i, yiy2i;
        size_t i = 0;
        while (xi >= x_min && xi + stepi < x_max && i <= N_max) {
            //LOG_DEBUG_CLI("Start RK4_SOE with localstecpcontrol", cfg);
            auto tmp1 = RK3_SOE_STEP(std::move(rhs1), std::move(rhs2), stepi, xi, vi, yi);
            auto tmp2 = RK3_SOE_STEP(std::move(rhs1), std::move(rhs2), stepi / 2, xi, vi, yi);
            tmp2 = RK3_SOE_STEP(std::move(rhs1), std::move(rhs2), stepi / 2, xi, std::get<0>(tmp2), std::get<1>(tmp2));
            Old_LocalError = LocalError;
            LocalError = std::max(std::abs(std::get<0>(tmp1) - std::get<0>(tmp2)), std::abs(std::get<1>(tmp1) - std::get<1>(tmp2))) / 3;
            //LOG_DEBUG_CLI("tmp generated", cfg);
            if (LocalError < eps / 16) {
                //LOG_DEBUG_CLI("try LocalError < lec", cfg);
                yi = std::get<1>(tmp1);
                vi = std::get<0>(tmp1);
                xi = xi + stepi;
                viv2i = std::abs(vi - std::get<0>(tmp2));
                yiy2i = std::abs(yi - std::get<1>(tmp2));

                i++;
                table.push_back({ xi, vi, yi, std::get<0>(tmp1), std::get<1>(tmp2), viv2i, yiy2i, LocalError, LocalError / Old_LocalError, stepi, C1, C2, 0.f, 0.f });
                //LOG_DEBUG_CLI("push_back done", cfg);
                stepi = stepi * 2;
                C2++;
                //LOG_DEBUG_CLI("LocalError < lec done", cfg);
            } else if (LocalError >= eps / 16 && LocalError <= eps) {
                //LOG_DEBUG_CLI("LocalError >= lec", cfg);
                yi = std::get<1>(tmp1);
                vi = std::get<0>(tmp1);
                xi = xi + stepi;
                viv2i = std::abs(vi - std::get<0>(tmp2));
                yiy2i = std::abs(yi - std::get<1>(tmp2));

                i++;
                table.push_back({ xi, vi, yi, std::get<0>(tmp2), std::get<1>(tmp2), viv2i, yiy2i, LocalError * 4, LocalError / Old_LocalError, stepi, C1, C2, 0.f, 0.f });
            } else if (LocalError > eps) {
                LOG_DEBUG_CLI("LocalError > lec*", cfg);
                stepi = stepi / 2;
                C1++;
            }
        }
        return table;
    }
    return (table);
}