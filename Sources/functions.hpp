#ifndef FUNCTIONS_HPP
#define FUNCTIONS_HPP
#include <cmath>
#include <QMap>
#include <QPair>
#include <QVector>
#include <algorithm>
#include <iostream>
class Functions
{
public:

    Functions(QMap<QString, double> _values);

    QVector<double> solver();

    QString getLog();

    void getSpeed(QVector<double> _roots, double _precision);
private:
    double a_0(QMap<QString, double> _values);
    double a_1(QMap<QString, double> _values);
    double a_2(QMap<QString, double> _values);
    double a_3(QMap<QString, double> _values);
    double a_4(QMap<QString, double> _values);
    double a_5(QMap<QString, double> _values);
    double a_6(QMap<QString, double> _values);


    int changesOfSign(QVector<double> _koef);

    QVector<double> signSwap(QVector<double> _koef);

    //число перемен знаков в ряде функции и ее производных
    int BF_series(QVector<double> _koef, double _point);

    //рассчет функции  в точке
    double funcResult(QVector<double> _koef, double _point);

    //n-я производная ряда
    QVector<double> diff(QVector<double> _koef, int _degree);

    void logWrite(QString _text, double _var);

    QVector<QPair<double, double>> splitSegment(QPair<double, double> _segment, int _parts);

    QVector<QPair<double, double>> getSegments(QPair<double, double> _segment, int depth);


public:
    double eps = 0.00001;
private:
    QVector<double> koef;

    QString log;

    QMap<QString, double> values;
};

#endif // FUNCTIONS_HPP
