#include <bits/stdc++.h>
using namespace std;

double f(double x, vector<double> &v) {
    double result = 0;
    int n = v.size();
    for (int i = 0; i < n; i++) {
        int deg = n - i - 1;
        result += v[i] * pow(x, deg);
    }
    return result;
}

// False Position for one bracket
pair<double, int> falsePosition(double a, double b, vector<double> &eqn, double tol, int maxIt = 1000) {
    double fa = f(a, eqn);
    double fb = f(b, eqn);
    if (fa * fb >= 0) return {0, 0};

    double c = a;
    int steps = 0;
    while (steps < maxIt) {
        c = b - fb * (b - a) / (fb - fa);
        double fc = f(c, eqn);

        steps++; // increment iteration

        if (fabs(fc) < tol) break;

        if (fa * fc < 0) {
            b = c;
            fb = fc;
        } else {
            a = c;
            fa = fc;
        }
    }
    return {c, steps};
}

int main() {
    ifstream in("F_input.txt");
    ofstream out("F_output.txt");

    if (!in || !out) return 1;

    int n;
    in >> n;
    double stepSize, tol;
    in >> stepSize >> tol;

    vector<double> eqn(n + 1);
    for (int i = 0; i <= n; i++) in >> eqn[i];

    out << "False Position Method\n";
    out << "Equation: ";
    for (int i = 0; i <= n; i++) {
        if (eqn[i] == 0) continue;
        if (i != 0 && eqn[i] > 0) out << " + ";
        if (eqn[i] < 0) out << " - ";
        out << abs(eqn[i]);
        int deg = n - i;
        if (deg > 0) out << "x";
        if (deg > 1) out << "^" << deg;
    }
    out << " = 0\n\n";

    // Fixed scanning range
    double scanStart = -10, scanEnd = 10;
    vector<pair<double, double>> intervals;

    double prev = scanStart;
    for (double cur = scanStart + stepSize; cur <= scanEnd; cur += stepSize) {
        if (f(prev, eqn) * f(cur, eqn) < 0)
            intervals.push_back({prev, cur});
        prev = cur;
    }

    if (intervals.empty()) {
        out << "No brackets found in this range\n";
        return 0;
    }

    out << left << setw(10) << "Root" 
        << setw(12) << "Iterations" 
        << setw(20) << "Bracket" << endl;
    out << "----------------------------------------\n";

    for (int i = 0; i < intervals.size(); i++) {
        auto res = falsePosition(intervals[i].first, intervals[i].second, eqn, tol);
        out << setw(10) << fixed << setprecision(6) << res.first
            << setw(12) << res.second
            << "[" << intervals[i].first << ", " << intervals[i].second << "]" << endl;
    }

    cout << "Output written to FB_output.txt\n";
    return 0;
}
