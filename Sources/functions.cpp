#include "functions.hpp"

QVector<QPair<double, double>> Functions::splitSegment(QPair<double, double> _segment, int _parts)
{
    QVector<QPair<double, double>> segments;
    double freq = (_segment.second - _segment.first) / _parts;
    for(int i = 0; i < _parts; i++)
    {
        segments.push_back(qMakePair<double, double>(_segment.first + i * freq, _segment.first + (i+1) * freq));
    }
    return segments;
}

QVector<QPair<double, double>> Functions::getSegments(QPair<double, double> _segment, int depth = 0)
{
    QVector<QPair<double, double>> temp;
    if(depth >= 8)
    {
        QVector<QPair<double, double>> temp;
        return temp;
    }
    else
    if((BF_series(koef, _segment.first) - BF_series(koef, _segment.second)) > 1)
    {
        double freq = double(_segment.second - _segment.first) / 10.0;
        for(int i = 0; i < 10; i++)
        {
            QVector<QPair<double, double>> rec = getSegments(qMakePair<double, double>(_segment.first + i * freq, _segment.first + (i+1) * freq), depth + 1);
            temp.append(rec);
        }
    }
    else
    if((BF_series(koef, _segment.first) - BF_series(koef, _segment.second)) == 0)
    {
        return QVector<QPair<double, double>>();
    }
    else
    if((BF_series(koef, _segment.first) - BF_series(koef, _segment.second)) == 1)
    {
        QVector<QPair<double, double>> result;
        result.push_back(_segment);
        return result;
    }
    if((BF_series(koef, _segment.first) - BF_series(koef, _segment.second)) < 0)
    {
        return QVector<QPair<double, double>>();
    }
    return temp;
}

void Functions::logWrite(QString _text, double _var)
{
    log.append(_text);
    log.append(QString::number(_var));
    log.append("\n");
}

int Functions::BF_series(QVector<double> _koef, double _point)
{
    QVector<double> temp;
    int max_degree = _koef.length();
    for(int i = 0; i < max_degree; ++i)
    {
        temp.push_back(funcResult(diff(_koef, i), _point));
    }
    int amount = changesOfSign(temp);
    return amount;
}

QVector<double> Functions::diff(QVector<double> _koef, int _degree)
{
    if(_degree == 0)
    {
        return _koef;
    }
    else if(_degree == 1)
    {
        QVector<double> temp;
        int max_degree = _koef.length() - 1;
        for(int i = 0; i < max_degree; ++i)
        {
            double temp_koef = (max_degree - i) * _koef.at(i);
            temp.push_back(temp_koef);
        }
        return temp;
    }
    else
    {
        QVector<double> temp;
        int max_degree = _koef.length() - 1;
        for(int i = 0; i < max_degree; ++i)
        {
            double temp_koef = (max_degree - i) * _koef.at(i);
            temp.push_back(temp_koef);
        }
        return diff(temp, _degree - 1);
    }

}

double Functions::funcResult(QVector<double> _koef, double _point)
{
    double result = 0;
    int max_degree = _koef.length() - 1;
    for(int i = 0; i < _koef.length(); ++i)
    {
        result = result +  _koef.at(i) * pow(_point, max_degree - i);
    }
    return result;
}

QString Functions::getLog()
{
    return log;
}

int Functions::changesOfSign(QVector<double> _koef)
{
    int amount = 0;
    for(int i = 1; i < _koef.length(); ++i)
    {
        if((_koef.at(i) * _koef.at(i - 1)) < 0.0)
        {
            amount ++;
        }
    }
    return amount;
}

QVector<double> Functions::signSwap(QVector<double> _koef)
{
    int max_degree = _koef.length() - 1;
    QVector<double> temp;

    for(int i = 0; i < _koef.length(); ++i)
    {
        if(((max_degree - i)%2) == 0)
        {
            temp.push_back(_koef.at(i));
        }
        else
        {
            temp.push_back(-_koef.at(i));
        }
    }
    return temp;
}

Functions::Functions(QMap<QString, double> _values)
{
    //инициализация констант уравнения
    koef.push_back(this->a_0(_values));
    log.append("a_0 был посчитан: ");
    log.append(QString::number(koef.back()));
    log.append("\n");

    koef.push_back(this->a_1(_values));
    log.append("a_1 был посчитан: ");
    log.append(QString::number(koef.back()));
    log.append("\n");

    koef.push_back(this->a_2(_values));
    log.append("a_2 был посчитан: ");
    log.append(QString::number(koef.back()));
    log.append("\n");

    koef.push_back(this->a_3(_values));
    log.append("a_3 был посчитан: ");
    log.append(QString::number(koef.back()));
    log.append("\n");

    koef.push_back(this->a_4(_values));
    log.append("a_4 был посчитан: ");
    log.append(QString::number(koef.back()));
    log.append("\n");

    koef.push_back(this->a_5(_values));
    log.append("a_5 был посчитан: ");
    log.append(QString::number(koef.back()));
    log.append("\n");

    koef.push_back(this->a_6(_values));
    log.append("a_6 был посчитан: ");
    log.append(QString::number(koef.back()));
    log.append("\n");

    values = _values;
}

double Functions::a_0(QMap<QString, double> _values)
{
    //константы необходимые для вычисления коэффициента
    double alpha_0 = (_values["gamma_0"] + 1.0) / (_values["gamma_0"] - 1.0);
    double alpha_3 = (_values["gamma_3"] + 1.0) / (_values["gamma_3"] - 1.0);
    double x = _values["p_3"] / _values["p_0"];
    double c_0 = sqrt(_values["gamma_0"] * _values["p_0"] / _values["rho_0"]);
    double c_3 = sqrt(_values["gamma_3"] * _values["p_3"] / _values["rho_3"]);
    double e_0 = (2 * pow(c_0, 2)) / (_values["gamma_0"] * (_values["gamma_0"] - 1) * pow((_values["u_3"] - _values["u_0"]), 2));
    double e_3 = (2 * pow(c_3, 2)) / (_values["gamma_3"] * (_values["gamma_3"] - 1) * pow((_values["u_3"] - _values["u_0"]), 2));

    //расчет коэффициента
    double a = pow((alpha_0 * e_3 - alpha_3 * x * e_0), 2);
    //добавлене коэффициента в список коэффициентов
    return a;
}

double Functions::a_1(QMap<QString, double> _values)
{
    //константы необходимые для вычисления коэффициента
    double alpha_0 = (_values["gamma_0"] + 1.0) / (_values["gamma_0"] - 1.0);
    double alpha_3 = (_values["gamma_3"] + 1.0) / (_values["gamma_3"] - 1.0);
    double x = _values["p_3"] / _values["p_0"];
    double c_0 = sqrt(_values["gamma_0"] * _values["p_0"] / _values["rho_0"]);
    double c_3 = sqrt(_values["gamma_3"] * _values["p_3"] / _values["rho_3"]);
    double e_0 = (2 * pow(c_0, 2)) / (_values["gamma_0"] * (_values["gamma_0"] - 1) * pow((_values["u_3"] - _values["u_0"]), 2));
    double e_3 = (2 * pow(c_3, 2)) / (_values["gamma_3"] * (_values["gamma_3"] - 1) * pow((_values["u_3"] - _values["u_0"]), 2));

    //рассчет коэффициента
    double a = 2 * ((alpha_0 * e_3 - alpha_3 * x * e_0) * (e_3 * (1 - 2 * alpha_0 * x) - e_0 * x * (x - 2 * alpha_3)) -
                         alpha_0 * alpha_3 * x * (alpha_0 * e_3 + alpha_3 * x * e_0));
    return a;
}

double Functions::a_2(QMap<QString, double> _values)
{
    //константы необходимые для вычисления коэффициента
    double alpha_0 = (_values["gamma_0"] + 1.0) / (_values["gamma_0"] - 1.0);
    double alpha_3 = (_values["gamma_3"] + 1.0) / (_values["gamma_3"] - 1.0);
    double x = _values["p_3"] / _values["p_0"];
    double c_0 = sqrt(_values["gamma_0"] * _values["p_0"] / _values["rho_0"]);
    double c_3 = sqrt(_values["gamma_3"] * _values["p_3"] / _values["rho_3"]);
    double e_0 = (2 * pow(c_0, 2)) / (_values["gamma_0"] * (_values["gamma_0"] - 1) * pow((_values["u_3"] - _values["u_0"]), 2));
    double e_3 = (2 * pow(c_3, 2)) / (_values["gamma_3"] * (_values["gamma_3"] - 1) * pow((_values["u_3"] - _values["u_0"]), 2));

    //рассчет коэффициента
    double a = pow(e_3, 2) * (6 * pow(alpha_0 * x, 2) - 8 * alpha_0 * x + 1) -
            2 * e_0 * e_3 * x * (alpha_0 * alpha_3 * (x*x + 4 * x + 1) - 2 * (x+1) * (alpha_3 + alpha_0 * x) + x) +
            pow(e_0 * x, 2) * (6 * alpha_3*alpha_3 - 8 * alpha_3 *x + x*x) + pow(alpha_0 * alpha_3 * x, 2) -
            2 * alpha_0 * x * e_3 * (alpha_0 * x - 2 * alpha_0 * alpha_3 * x + 2 * alpha_3) -
            2 * alpha_3 * x * x * e_0 * (alpha_3 + 2 * alpha_0 * x - 2 * alpha_0 * alpha_3);

    return a;
}

double Functions::a_3(QMap<QString, double> _values)
{
    //константы необходимые для вычисления коэффициента
    double alpha_0 = (_values["gamma_0"] + 1.0) / (_values["gamma_0"] - 1.0);
    double alpha_3 = (_values["gamma_3"] + 1.0) / (_values["gamma_3"] - 1.0);
    double x = _values["p_3"] / _values["p_0"];
    double c_0 = sqrt(_values["gamma_0"] * _values["p_0"] / _values["rho_0"]);
    double c_3 = sqrt(_values["gamma_3"] * _values["p_3"] / _values["rho_3"]);
    double e_0 = (2 * pow(c_0, 2)) / (_values["gamma_0"] * (_values["gamma_0"] - 1) * pow((_values["u_3"] - _values["u_0"]), 2));
    double e_3 = (2 * pow(c_3, 2)) / (_values["gamma_3"] * (_values["gamma_3"] - 1) * pow((_values["u_3"] - _values["u_0"]), 2));

    //рассчет коэффициента
    double a = -2 * x * (2 * e_3 * e_3 * (pow(alpha_0 * x, 2) - 3 * alpha_0 * x + 1) + e_0 * e_3 * ((alpha_3 + alpha_0 * x) *
            (x*x + 4 * x + 1) - 2 * alpha_0 * alpha_3 * x * (x + 1) - 2 * x * (x + 1)) +
            2 * e_0*e_0 * x * (x*x - 3 * alpha_3 * x + alpha_3*alpha_3) - alpha_0 * alpha_3 * x * (alpha_0 * x + alpha_3) +
            e_3 * (alpha_0*alpha_0 * alpha_3 * x*x - 2 * x * (2 * alpha_0 * alpha_3 + alpha_0*alpha_0 * x) + (2 * alpha_0 * x + alpha_3)) +
            e_0 * x * (alpha_0 * alpha_3*alpha_3 - 2 * alpha_3 * (alpha_3 + 2 * alpha_0 * x) + 2 * alpha_3 * x + alpha_0 * x*x));

    return a;
}


double Functions::a_4(QMap<QString, double> _values)
{
    //константы необходимые для вычисления коэффициента
    double alpha_0 = (_values["gamma_0"] + 1.0) / (_values["gamma_0"] - 1.0);
    double alpha_3 = (_values["gamma_3"] + 1.0) / (_values["gamma_3"] - 1.0);
    double x = _values["p_3"] / _values["p_0"];
    double c_0 = sqrt(_values["gamma_0"] * _values["p_0"] / _values["rho_0"]);
    double c_3 = sqrt(_values["gamma_3"] * _values["p_3"] / _values["rho_3"]);
    double e_0 = (2 * pow(c_0, 2)) / (_values["gamma_0"] * (_values["gamma_0"] - 1) * pow((_values["u_3"] - _values["u_0"]), 2));
    double e_3 = (2 * pow(c_3, 2)) / (_values["gamma_3"] * (_values["gamma_3"] - 1) * pow((_values["u_3"] - _values["u_0"]), 2));

    //рассчет коэффициента
    double a = x*x * ( e_3*e_3 * (pow(alpha_0 * x, 2) - 8 * alpha_0 * x + 6) -
            2 * e_0 * e_3 * (alpha_0 * alpha_3 * x - 2 * (x+1) * (alpha_3 + alpha_0 * x) + x*x + 4 * x + 1) +
            e_0*e_0 * (alpha_3*alpha_3 - 8 * alpha_3 * x + 6 * x*x) + (alpha_3*alpha_3 + 4 * alpha_0 * alpha_3 * x + pow(alpha_0 * x, 2)) -
            2 * e_3 * ((alpha_0*alpha_0 * x + 2 * alpha_0 * alpha_3) * x - 2 * (2 * alpha_0 * x + alpha_3) + 1) -
            2 * e_0 * (alpha_3 * (2 * alpha_0 * x + alpha_3) - 2 * x * (2 * alpha_3 + alpha_0 * x) + x*x));

    return a;
}

double Functions::a_5(QMap<QString, double> _values)
{
    //константы необходимые для вычисления коэффициента
    double alpha_0 = (_values["gamma_0"] + 1.0) / (_values["gamma_0"] - 1.0);
    double alpha_3 = (_values["gamma_3"] + 1.0) / (_values["gamma_3"] - 1.0);
    double x = _values["p_3"] / _values["p_0"];
    double c_0 = sqrt(_values["gamma_0"] * _values["p_0"] / _values["rho_0"]);
    double c_3 = sqrt(_values["gamma_3"] * _values["p_3"] / _values["rho_3"]);
    double e_0 = (2 * pow(c_0, 2)) / (_values["gamma_0"] * (_values["gamma_0"] - 1) * pow((_values["u_3"] - _values["u_0"]), 2));
    double e_3 = (2 * pow(c_3, 2)) / (_values["gamma_3"] * (_values["gamma_3"] - 1) * pow((_values["u_3"] - _values["u_0"]), 2));

    //рассчет коэффициента
    double a = 2 * pow(x, 3) * (e_3*e_3 * (alpha_0 * x - 2) - e_0 * e_3 * (alpha_0 * x - 2 + alpha_3 - 2 * x) +
            e_0*e_0 * (alpha_3 - 2 * x) + (alpha_3 + alpha_0 * x) -
            e_3 * (2 * alpha_0 * x + alpha_3 - 2) - e_0 * (2 * alpha_3 + alpha_0 * x - 2 * x));

    return a;
}

double Functions::a_6(QMap<QString, double> _values)
{
    //константы необходимые для вычисления коэффициента
    double alpha_0 = (_values["gamma_0"] + 1.0) / (_values["gamma_0"] - 1.0);
    double alpha_3 = (_values["gamma_3"] + 1.0) / (_values["gamma_3"] - 1.0);
    double x = _values["p_3"] / _values["p_0"];
    double c_0 = sqrt(_values["gamma_0"] * _values["p_0"] / _values["rho_0"]);
    double c_3 = sqrt(_values["gamma_3"] * _values["p_3"] / _values["rho_3"]);
    double e_0 = (2 * pow(c_0, 2)) / (_values["gamma_0"] * (_values["gamma_0"] - 1) * pow((_values["u_3"] - _values["u_0"]), 2));
    double e_3 = (2 * pow(c_3, 2)) / (_values["gamma_3"] * (_values["gamma_3"] - 1) * pow((_values["u_3"] - _values["u_0"]), 2));

    //рассчет коэффициента
    double a = pow(x, 4) * (pow(e_3 - e_0, 2) + 1 - 2 * (e_3 + e_0));

    return a;
}

QVector<double> Functions::solver()
{
    //локализация корней
    log.append("Локализация корней \n");
    QVector<double> tempA(koef);
    QVector<double> tempB(koef);
    tempA.pop_front();
    tempB.pop_back();
    double A = abs(tempA.at(0));
    double B = abs(tempB.at(0));
    for(int i = 0; i < (koef.length() - 1); ++i)
    {
        if(abs(tempA.at(i)) > A)
        {
            A = abs(tempA.at(i));
        }

        if(abs(tempB.at(i) > B))
        {
            B = abs(tempB.at(i));
        }
    }

    if(((abs(koef.back()) + B) == 0)||(koef.front() == 0))
    {
        log.append("Уравнение не корректно \n");
        return QVector<double>();

    }

    double z_min = abs(koef.back()) / (koef.back() + B);
    double z_max = 1 + A / abs(koef.front());

    log.append("z min: ");
    log.append(QString::number(z_min));
    log.append("\n");

    log.append("z max: ");
    log.append(QString::number(z_max));
    log.append("\n");

    //число положительных корней
    double pos = this->changesOfSign(koef);
    //число отрицательных корней
    double neg = this->changesOfSign(this->signSwap(koef));
    logWrite("Максимальное число положительных корней: ", pos);
    logWrite("Максимальное число отрицательных корней: ", neg);

    //Используем теорему Бюдана-Фурье
    log.append("Уточним корни с помощью теоремы Бюдана-Фурье\n");

    log.append("Разбиваем отрезки, на котрых могут находиться корни, на множество меньших отрезков, на концах которых, "
               "считаем число перемен знаков в ряде функции и ее производных \n");

    QVector<QPair<double, double>> segments;
    segments.append(getSegments(qMakePair<double, double>(-z_max, -z_min), 0));
    segments.append(getSegments(qMakePair<double, double>(z_min, z_max), 0));

    QVector<double> roots;
    //уточнение корней
    while(segments.length() > 0)
    {
        int depth = 0;

        QPair<double, double> seg = segments.takeFirst();
        if((funcResult(koef, seg.first) * funcResult(koef, seg.second)) > 0)
        {
           continue;
        }
        else
        {
            while(((seg.second - seg.first) > eps) && (depth < 10))
            {

                double midle = (seg.first + seg.second) / 2.0;
                if((funcResult(koef, seg.first) * funcResult(koef, midle)) < 0)
                {
                    seg.second = midle;
                }
                else if((funcResult(koef, seg.second) * funcResult(koef, midle)) < 0)
                {
                    seg.first = midle;
                }
                ++depth;
            }
        }
        if(depth < 100)
        {
            roots.push_back(seg.first);
        }

    }
    for(double i: roots)
    {
        logWrite("Найден корень:  ", i);
    }

    return roots;
}

void Functions::getSpeed(QVector<double> roots, double precision)
{
    QVector<double> d0;
    QVector<double> d3;
    QVector<double> _roots;
    for(int t = 0; t < roots.length(); ++t)
    {
        if(!(roots.at(t) < 0))
        {
            _roots.push_back(roots.at(t)* values["p_0"]);
        }
    }

    QVector<double> u;
    QVector<double> p_checked;

    for(int t = 0; t < _roots.length(); ++t)
    {
        double c_0 = sqrt(values["gamma_0"] * values["p_0"] / values["rho_0"]);
        double c_3 = sqrt(values["gamma_3"] * values["p_3"] / values["rho_3"]);
        // U_1 options

        double tempU1_neg = values["u_0"] - ((_roots.at(t) - values["p_0"]) / (values["rho_0"] * c_0)) * sqrt((2 * values["gamma_0"]) / ((values["gamma_0"] - 1) + ((values["gamma_0"] + 1) * _roots.at(t) / values["p_0"])));
        double tempU1_pos = values["u_0"] + ((_roots.at(t) - values["p_0"]) / (values["rho_0"] * c_0)) * sqrt((2 * values["gamma_0"]) / ((values["gamma_0"] - 1) + ((values["gamma_0"] + 1) * _roots.at(t) / values["p_0"])));
        // U_2 options
        double tempU2_neg = values["u_3"] - ((_roots.at(t) - values["p_3"]) / (values["rho_3"] * c_3)) * sqrt((2 * values["gamma_3"]) / ((values["gamma_3"] - 1) + ((values["gamma_3"] + 1) * _roots.at(t) / values["p_3"])));
        double tempU2_pos = values["u_3"] + ((_roots.at(t) - values["p_3"]) / (values["rho_3"] * c_3)) * sqrt((2 * values["gamma_3"]) / ((values["gamma_3"] - 1) + ((values["gamma_3"] + 1) * _roots.at(t) / values["p_3"])));

        // looking for U_1 = U_2 (from the relation at the contact discontinuity) ! with the same P
        if ((abs(tempU1_neg - tempU2_neg) < precision) || (abs(tempU1_neg - tempU2_pos) < precision))
        {
            u.push_back(tempU1_neg);
            p_checked.push_back(_roots.at(t));
            logWrite("Найдена скорость U_1 ", tempU1_neg);
        }
        if ((abs(tempU1_pos - tempU2_neg) < precision) || (abs(tempU1_pos - tempU2_pos) < precision))
        {
            u.push_back(tempU1_pos);
            p_checked.push_back(_roots.at(t));
            logWrite("Найдена скорость U_1 ", tempU1_pos);
        }
    }
    for(int i = 0; i < p_checked.length(); ++i)
    {
        d0.push_back(u.at(i) + (p_checked.at(i) - values["p_0"]) / (values["rho_0"] * (u.at(i) - values["u_0"])));
        d3.push_back(u.at(i) + (p_checked.at(i) - values["p_3"]) / (values["rho_3"] * (u.at(i) - values["u_3"])));
        logWrite("Найдена скорость D_0 ", d0.last());
        logWrite("Найдена скорость D_3 ", d3.last());
    }

}
