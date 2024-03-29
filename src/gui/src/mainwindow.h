#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include "helpform.h"
#include <QWidget>
#include <QVector>
#include <cmath>
#include <cstddef>

#include "nmlib.hpp"

QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();
    

private slots:
    void on_button_plot_clicked();

    void on_button_clear_clicked();

    void on_button_save_points_clicked();

    void on_button_plot_from_file_clicked();

    void on_exit_button_clicked();

    void on_getdata_buttom_clicked();

    void on_radioButton_blue_clicked(bool checked);

    void on_radioButton_red_clicked(bool checked);

    void on_radioButton_green_clicked(bool checked);

    void on_radioButton_violet_clicked(bool checked);

    void on_radioButton_mistake_clicked(bool checked);

    void on_button_table_clicked();

    void on_HelpButton_clicked();

private:
    Ui::MainWindow *ui;
    double h, X;
    bool LEC = false;
    double x_begin, x_end;
    double x_start, y_start;
    double precision;
    uint N;
    int count_plot = 0;
    double du;
    resultTable res1;
    double A, B, C;   
    double brd_prec = 0.l;
    QColor col = QColor(0, 0, 255);

    // double max_LE = 0.l;
    // double max_LE_x = 0.l;
    // double max_step = 0.l;
    // double max_step_x = 0.l;
    // double min_step = 0.l;
    // double min_step_x = 0.l;
    // double max_uvi = 0.l;
    // double max_uvi_x = 0.l;
    // double last_x = 0.l;
    // double last_v = 0.l;
    // size_t steps_num = 0;
    // size_t total_step_augs = 0;
    // size_t total_step_dims = 0;
    config cfg;

    HelpForm form;
};
#endif // MAINWINDOW_H
