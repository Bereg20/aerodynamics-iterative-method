#include <iostream>
#include <cmath>
#include <fstream>
#define pi 3.14159265358979323846
using namespace std;

struct Constant {
    double f = 200;
    double omega = 2 * pi * f;
    double l = 0.5 / 100.0;
    double L = 6.5 / 100.0;
    double a = 0.5 / 100.0;
    double ro = 1.117;
    double F = 2 * pi * pow((L - l), 2);

    double alpha_0 = 3 * pi / 180.0;
};

// Нулевая итерация (без W)
double LiftZero(double alpha, Constant& params) {
    int N = 1000;
    double total_force = 0.0;
    double h = (params.L - params.l) / N;

    for (int i = 0; i < N; i++) {
        double r = params.l + i * h + h / 2.0;
        double V_local = params.omega * r;
        double dL = 2.0 * alpha * params.ro * V_local * V_local * pi * params.a * h;
        total_force += dL;
    }
    return total_force;
}

// Итерационная функция (с учётом W)
double LiftWithW(double alpha_install, double W_induced, Constant& params) {
    int N = 1000;
    double total_force = 0.0;
    double h = (params.L - params.l) / N;

    for (int i = 0; i < N; i++) {
        double r = params.l + i * h;
        double V_okr = params.omega * r;
        double V_full = sqrt(V_okr * V_okr + W_induced * W_induced);
        double beta_local = atan(W_induced / V_okr);
        double alpha_eff = alpha_install - beta_local;
        double dL = 2.0 * alpha_eff * params.ro * V_full * V_full * pi * params.a * h;
        total_force += dL;
    }
    return total_force;
}

double W(double P, Constant& params) {
    return sqrt(P / (2 * params.F * params.ro));
}

int main() {
    ofstream file;
    file.open("iteration.txt");

    setlocale(LC_ALL, "Russian");
    Constant c;
    const double epsilon = 1e-6;

    // Нулевая итерация
    double alpha_install = c.alpha_0;
    double Force = LiftZero(alpha_install, c);
    double V_podsos = W(Force, c);

    cout << "Исходные данные: " << endl;
    cout << "Угол установки: " << alpha_install << " рад (" << alpha_install * 180 / pi << "°)" << endl;
    cout << "Тяга (0 итерация): " << Force << " Н" << endl;
    cout << "Индуцированная скорость: " << V_podsos << " м/с" << endl;
    cout << "\nВведите количество итераций: ";
    // cin >> n;

    for (int i = 1; i <= 1000; i++) {
        double force_old = Force;

        Force = LiftWithW(alpha_install, V_podsos, c);
        V_podsos = W(Force, c);

        cout << "Итерация " << i << ": Force = " << Force << " Н, W = " << V_podsos << " м/с" << endl;
        file << i << "    " << Force << endl;

        if (abs(force_old - Force) < epsilon) {
            cout << "Процесс сошелся" << endl;
            break;
        }
    }

    file.close();
    cout << "\nРезультат: Force = " << Force << " Н" << endl;

    return 0;
}