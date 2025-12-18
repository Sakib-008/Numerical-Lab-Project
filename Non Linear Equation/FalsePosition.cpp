#include <bits/stdc++.h>
using namespace std;

double f(double x, const vector<double>& a) {
    int n = a.size() - 1;
    double res = 0;
    for (int i = 0; i <= n; i++)
        res += a[i] * pow(x, n - i);
    return res;
}
// Override the func to change function
pair<double,int> falsePos(double l, double r, const vector<double>& a,
                          double eps, int maxIt = 1000) {
    if (f(l,a) * f(r,a) >= 0) return {0.0, 0};

    double x = l;
    int it = 0;

    while (it < maxIt && fabs(r - l) > eps) {
        it++;
        double fl = f(l,a), fr = f(r,a);
        x = r - fr * (r - l) / (fr - fl);
        double fx = f(x,a);

        if (fl * fx < 0) r = x;
        else l = x;
    }
    return {x, it};
}

vector<pair<double,double>> findBrackets(const vector<double>& a,
                                         double xmin, double xmax, double h) {
    vector<pair<double,double>> b;
    for (double x = xmin; x < xmax; x += h) {
        if (f(x,a) * f(x + h,a) < 0)
            b.push_back({x, x + h});
    }
    return b;
}

int main() {
    ifstream in("FalsePositionIn.txt");
    ofstream out("FalsePositionOut.txt");
    if (!in || !out) return 1;

    int n;
    in >> n;
    vector<double> a(n);
    for (int i = 0; i < n; i++) in >> a[i];

    double step, eps;
    in >> step >> eps;

    double xmax = sqrt(pow(a[1]/a[0], 2) - 2*(a[2]/a[0]));
    double xmin = -xmax;

    auto brackets = findBrackets(a, xmin, xmax, step);

    out << fixed << setprecision(4);
    out << "Root      Iter   Bracket\n";
    out << "------------------------------\n";
    cout<<"Output in FalsePositionOut.txt";

    for (auto &br : brackets) {
        auto res = falsePos(br.first, br.second, a, eps);
        out << setw(8) << res.first
            << setw(8) << res.second
            << "   [" << br.first << ", " << br.second << "]\n";
    }

    return 0;
}
