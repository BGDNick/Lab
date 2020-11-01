#include "mainwindow.hpp"
#include "ui_mainwindow.h"

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    //generating values from task
    gamma_0 = 5.0 / 3.0;
    gamma_3 = 5.0 / 3.0;//
    rho_0 = 7.9;
    rho_1 = 0;
    rho_2 = 0;
    rho_3 = 11.37;
    u_0 = 0;
    u_3 = 5 * pow(10, 4);
    p_0 = 3.04 * pow(10, 9);
    p_3 = 1.17928 *  pow(10, 9);
    d_0 = 0;
    d_3 = 0;

    values["gamma_0"] = gamma_0;
    values["gamma_3"] = gamma_3;
    values["rho_0"] = rho_0;
    values["rho_1"] = rho_1;
    values["rho_2"] = rho_2;
    values["rho_3"] = rho_3;
    values["u_0"] = u_0;
    values["u_3"] = u_3;
    values["p_0"] = p_0;
    values["p_3"] = p_3;
    values["d_0"] = d_0;
    values["d_3"] = d_3;

    ui->doubleSpinBoxGamma_0->setValue(gamma_0);
    ui->doubleSpinBoxGamma_3->setValue(gamma_3);
    ui->doubleSpinBoxP_0->setValue(p_0);
    ui->doubleSpinBoxP_3->setValue(p_3);
    ui->doubleSpinBoxRho_0->setValue(rho_0);
    ui->doubleSpinBoxRho_3->setValue(rho_3);
    ui->doubleSpinBoxU_0->setValue(u_0);
    ui->doubleSpinBoxU_3->setValue(u_3);
    this->setFixedSize(800, 400);
}

MainWindow::~MainWindow()
{
    delete ui;
}


void MainWindow::on_doubleSpinBoxGamma_0_valueChanged(double arg1)
{
    gamma_0 = arg1;
}

void MainWindow::on_doubleSpinBoxRho_0_valueChanged(double arg1)
{
    rho_0 = arg1;
}

void MainWindow::on_doubleSpinBoxU_0_valueChanged(double arg1)
{
    u_0 = arg1;
}

void MainWindow::on_doubleSpinBoxP_0_valueChanged(double arg1)
{
    p_0 = arg1;
}

void MainWindow::on_doubleSpinBoxGamma_3_valueChanged(double arg1)
{
    gamma_3 = arg1;
}

void MainWindow::on_doubleSpinBoxRho_3_valueChanged(double arg1)
{
    rho_3 = arg1;
}

void MainWindow::on_doubleSpinBoxU_3_valueChanged(double arg1)
{
    u_3 = arg1;
}

void MainWindow::on_doubleSpinBoxP_3_valueChanged(double arg1)
{
    p_3 = arg1;
}

void MainWindow::on_pushButton_clicked()
{

    values["gamma_0"] = gamma_0;
    values["gamma_3"] = gamma_3;
    values["rho_0"] = rho_0;
    values["rho_1"] = rho_1;
    values["rho_2"] = rho_2;
    values["rho_3"] = rho_3;
    values["u_0"] = u_0;
    values["u_3"] = u_3;
    values["p_0"] = p_0;
    values["p_3"] = p_3;
    values["d_0"] = d_0;
    values["d_3"] = d_3;

    QVector<double> koef = {0.971332, -9.07362, 17.3308, -3.37955, 0.0245529, 0.0137498, 0.00028348};
    Functions f(values);
    f.eps = ui->doubleSpinBoxEps->value();
    f.getSpeed(f.solver(), 100);
    ui->plainTextEdit->setPlainText(f.getLog());
}
