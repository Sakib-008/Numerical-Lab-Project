#include <bits/stdc++.h>
using namespace std;

double f(double x) {
    return x*x - 4.0;
}

double df(double x) {
    return 2.0*x;
}
double secant(double x_n, double x_n1, double epsilon, int max_iteration) {
    double x_n2;
    for (int i = 0; i < max_iteration; i++) {
        x_n2 = x_n1 - ((x_n1 - x_n) / (f(x_n1) - f(x_n))) * f(x_n1);
        if (abs(x_n2 - x_n1) < epsilon) {
            return x_n2;
        }
        x_n = x_n1;
        x_n1 = x_n2;
    }
    return NAN;
}

int main() {
    ifstream fin("input.txt");
    ofstream fout("output.txt");

    double range_start, range_end, step, epsilon;
    int max_iteration;
    fin >> range_start >> range_end >> step >> epsilon >> max_iteration;

    vector<double> raphson_points;
    vector<pair<double, double>> secant_pairs;

    for (double i = range_start; i < range_end; i += step) {
        if (f(i) * f(i + step) < 0) {
            raphson_points.push_back(i + step / 2);
            secant_pairs.push_back({i, i + step});
        }
    }
    if (raphson_points.size() == 0) {
        raphson_points.push_back(0.5);
    }

    for (auto pr : secant_pairs) {
        double root = secant(pr.first, pr.second, epsilon, max_iteration);
        fout << "Root using Secant: " << root << endl;
    }

    fin.close();
    fout.close();
    return 0;
}
