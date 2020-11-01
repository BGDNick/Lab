#ifndef MAINWINDOW_HPP
#define MAINWINDOW_HPP
#include <cmath>
#include "functions.hpp"
#include <iostream>
#include <QMap>
#include <QMainWindow>
#include <QVector>
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
    void on_doubleSpinBoxGamma_0_valueChanged(double arg1);

    void on_doubleSpinBoxRho_0_valueChanged(double arg1);

    void on_doubleSpinBoxU_0_valueChanged(double arg1);

    void on_doubleSpinBoxP_0_valueChanged(double arg1);

    void on_doubleSpinBoxGamma_3_valueChanged(double arg1);

    void on_doubleSpinBoxRho_3_valueChanged(double arg1);

    void on_doubleSpinBoxU_3_valueChanged(double arg1);

    void on_doubleSpinBoxP_3_valueChanged(double arg1);

    void on_pushButton_clicked();

private:
    Ui::MainWindow *ui;
private:
        double gamma_0;
        double gamma_3;
        double rho_0;
        double rho_1;
        double rho_2;
        double rho_3;
        double u_0;
        double u_3;
        double p_0;
        double p_3;
        double d_0;
        double d_3;

        QMap<QString, double> values;

};
#endif // MAINWINDOW_HPP
