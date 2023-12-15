#include "./ui_mainwindow.h"
#include "logger.hpp"
#include "mainwindow.h"
#include "nmlib.hpp"
#include <cmath>
#include <cstddef>
#include <fstream>
#include <qcustomplot.h>
#include <qvalidator.h>
#include <string>
#include <vector>

MainWindow::MainWindow(QWidget* parent) : QMainWindow(parent), ui(new Ui::MainWindow) {
    ui->setupUi(this);

    this->ui->plot->xAxis->setLabel("t");
    this->ui->plot->yAxis->setLabel("x(t)");

    this->ui->plot->xAxis->setRange(-10, 10);
    this->ui->plot->yAxis->setRange(-7, 7);

    this->ui->plot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);

    func = 0;
}

MainWindow::~MainWindow() {
    delete ui;
}

void MainWindow::on_button_plot_clicked() {
    QVector<double> x, v;

    x_begin = this->ui->lineEdit_x_start->text().toDouble();
    x_end = this->ui->lineEdit_x_end->text().toDouble();
    h = this->ui->lineEdit_step->text().toDouble();
    precision = this->ui->lineEdit_precision->text().toDouble();

    for (size_t i = 0; i < res1.size(); i++) {
        x.push_back(res1.at(i).xi);
        v.push_back(res1.at(i).vi);
    }

    this->ui->plot->addGraph();

    QPen pen;
    pen.setWidth(2);
    pen.setColor(col);
    this->ui->plot->graph(count_plot)->setPen(pen);
    this->ui->plot->addGraph();
    pen.setWidth(2);
    pen.setColor(col);

    this->ui->plot->graph(count_plot)->addData(x, v);
    this->ui->plot->replot();
    count_plot++;
}

void MainWindow::on_button_clear_clicked() {
    for (int i = 0; i < count_plot; i++) {
        this->ui->plot->graph(i)->data()->clear();
    }
    count_plot = 0;
    this->ui->plot->replot();
    this->ui->plot->update();
}

void MainWindow::on_button_save_points_clicked() {
    std::ofstream os("table.txt");

    os << "x\tu1\tu2\tu1_2\tu2_2\tu1-u1_2\tu2-u2_2\tОЛП\tОЛП/ОЛПП\th"
       << "C1\tC2\tu1-U1\tu2-U2\n";
    dumpTableInFile(res1, os);
}

void MainWindow::on_button_plot_from_file_clicked() {}

void MainWindow::on_exit_button_clicked() {
    QMessageBox::StandardButton exit = QMessageBox::question(this, " ", "you wanna exit?", QMessageBox::Yes | QMessageBox::No);
    if (exit == QMessageBox::Yes) {
        QApplication::quit();
    }
}

void MainWindow::on_getdata_buttom_clicked() {
    x_begin = this->ui->lineEdit_x_start->text().toDouble();
    x_end = this->ui->lineEdit_x_end->text().toDouble();
    h = this->ui->lineEdit_step->text().toDouble();
    precision = this->ui->lineEdit_precision->text().toDouble();
    x_start = this->ui->lineEdit_start_x->text().toDouble();
    y_start = this->ui->lineEdit_start_y->text().toDouble();
    du = this->ui->lineEdit_last->text().toDouble();
    N = this->ui->lineEdit_n->text().toInt();
    A = this->ui->lineEdit_a->text().toDouble();
    B = this->ui->lineEdit_b->text().toDouble();
    C = this->ui->lineEdit_c->text().toDouble();
    brd_prec = this->ui->lineEdit_bdr_pec->text().toDouble();

    cfg = { x_begin, x_end, x_start, y_start, du, h, N, LEC, precision, A, B, C, brd_prec, 0b11111111111111 };

    // for system w/ vi and yi w/o LEC (for LEC default value is applicable)
    if (!LEC) {
        cfg.set_col_visibility(0b11000000000011);
    } else {
        cfg.set_col_visibility(0b11111110000011);
    }
    // if you have only one variable u(x) in your task, you need to specify col_visibility as 0b11000000000011 w/ no lec
    // and 0b11111110101011 if lec is off (this both cases applicable only if you have a true solution. 
    // Otherwise you need to set two first bits to zero in cases discribed above)

    switch (func) {
    case 0:
        // res1 = task_rk4(test_rhs, cfg);
        break;
    case 1:
        // res1 = task_rk4(task1_rhs, cfg);
        break;
    case 2: {
        static auto task2_rhs1 = [&](double x, double u, double du) {
            return du;
        };
        static auto task2_ths2 = [&](double x, double u, double du) {
            return -(B / A) * du - (C / A) * u;
        };
        res1 = utils::RK3_SOE(std::move(task2_rhs1), std::move(task2_ths2), cfg);

        static const double bp = std::sqrt(4 * A * C - B * B) / (2 * A);
        static const double ap = -(B / (2 * A));
        static const double fp = -1 * ap / bp;
        calculate_global_error(res1, [&](const double& x) -> double {
            return 10 * std::exp(ap * x) * (std::cos(bp * x) + fp * std::sin(bp * x));
        });
        break;
    }
    default:
        break;
    }

    // max_LE = find_max_LE(res1);
    // max_step = find_max_h(res1);
    // max_uvi = find_max_uvi(res1);
    // min_step = find_min_h(res1);
    // steps_num = res1.size() - 1;

    auto [max_LE, max_LE_x] = find_max_LE(res1);
    auto [max_step, max_step_x] = find_max_h(res1);
    auto [min_step, min_step_x] = find_min_h(res1);
    auto [max_uvi, max_uvi_x] = find_max_uvi(res1);

    size_t steps_num = res1.size() - 1;
    size_t total_augs = res1.back().C2;
    size_t total_dims = res1.back().C1;

    double last_x = res1.back().xi;
    double last_v = res1.back().vi;

    double brd_dist = x_end -  res1.back().xi;

    if(!LEC) {
        this->ui->max_le_txt->setText("--");
        this->ui->max_ole_t->setText("--");
        this->ui->step_aug->setText("--");
        this->ui->step_dim->setText("--");
    } else {
        this->ui->max_le_txt->setText(QString::number(max_LE));
        this->ui->max_ole_t->setText(QString::number(max_LE_x));
        this->ui->step_aug->setText(QString::number(total_augs));
        this->ui->step_dim->setText(QString::number(total_dims));
    }

    this->ui->max_h_txt->setText(QString::number(max_step));
    this->ui->max_step_t->setText(QString::number(max_step_x));
    this->ui->min_h_txt->setText(QString::number(min_step));
    this->ui->min_step_t->setText(QString::number(min_step_x));
    this->ui->max_glob_txt->setText(QString::number(max_uvi));
    this->ui->max_glob_t->setText(QString::number(max_uvi_x));
    this->ui->step_num_txt->setText(QString::number(steps_num));
    this->ui->last_t->setText(QString::number(last_x));
    this->ui->last_v->setText(QString::number(last_v));
    this->ui->right_end_d->setText(QString::number(brd_dist));
}

//void MainWindow::on_Help_buttom_clicked(){};

void MainWindow::on_radioButton_blue_clicked(bool checked) {
    if (checked) {
        col = QColor(0, 0, 255);
    }
}

void MainWindow::on_radioButton_red_clicked(bool checked) {
    if (checked) {
        col = QColor(255, 0, 0);
    }
}

void MainWindow::on_radioButton_green_clicked(bool checked) {
    if (checked) {
        col = QColor(0, 255, 0);
    }
}

void MainWindow::on_radioButton_violet_clicked(bool checked) {
    if (checked) {
        col = QColor(105, 0, 198);
    }
}

void MainWindow::on_radioButton_mistake_clicked(bool checked) {
    LEC = checked;
}

void MainWindow::on_button_table_clicked() {
    ui->tableWidget->clear();

    size_t colsnum = 0;
    for (size_t i = 0; i < 14; i++) {
        colsnum = colsnum + static_cast<size_t>(cfg.isColVisible[i]);
    }
    ui->tableWidget->setRowCount(res1.size());
    ui->tableWidget->setColumnCount(colsnum);

    std::cout << res1.size() << std::endl;

    std::vector<const char*> horizontalLabels = { "t", "v", "y", "v_2", "v2_2", "v-v_2", "v2-v2_2", "ОЛП", "ОЛП/ОЛПП", "h", "Ум. шаг", "Ув. шаг", "u", "v-u" };
    QStringList hlabels;

    for (size_t i = 0; i < 14; i++) {
        if (cfg.is_col_visible(i)) {
            hlabels << horizontalLabels[i];
        }
    }
    ui->tableWidget->setHorizontalHeaderLabels(hlabels);

    for (size_t i = 0; i < res1.size(); i++) {
        auto row_tuple = res1.at(i).get_tuple();
        size_t j = 0, col = 0;
        apply_elemwise(
            [&](const auto& elem) {
                if (cfg.is_col_visible(j)) {
                    QTableWidgetItem* item = new QTableWidgetItem(QString::number(elem));
                    ui->tableWidget->setItem(i, col, item);
                    col++;
                }
                j++;
            },
            row_tuple,
            std::make_index_sequence<std::tuple_size<decltype(row_tuple)>::value>{});
    }
}

void MainWindow::on_comboBox_activated(int index) {
    func = index;
}

void MainWindow::on_HelpButton_clicked() {
    form.show();
}
