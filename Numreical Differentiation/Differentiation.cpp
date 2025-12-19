#include <bits/stdc++.h>
using namespace std;

double func(double x) {
    return x*x*x - 2*x*x + x;
}

void forwardDifferentiation(double a, double b, int n, double p, ofstream &fout) {
    double h = (b - a) / n;

    vector<double> X(n+1), Y(n+1);
    for (int i = 0; i <= n; i++) {
        X[i] = a + i*h;
        Y[i] = func(X[i]);
    }

    // Forward difference table
    vector<vector<double>> d(n+1, vector<double>(n+1, 0));
    for (int i = 0; i <= n; i++)
        d[i][0] = Y[i];

    for (int j = 1; j <= n; j++)
        for (int i = 0; i <= n - j; i++)
            d[i][j] = d[i+1][j-1] - d[i][j-1];

    fout << "Forward Difference Table:\n";
    for (int i = 0; i <= n; i++) {
        for (int j = 0; j <= n-i; j++)
            fout << d[i][j] << "\t";
        fout << "\n";
    }

    double u = (p - X[0]) / h;

    double fp = 0, fpp = 0;

    if (n >= 2) {
        fp = (d[0][1] + ((2*u - 1)/2.0)*d[0][2] + ((3*u*u - 6*u + 2)/6.0)*d[0][3]) / h;
        fpp = (d[0][2] + (u - 1)*d[0][3]) / (h*h);
    } else if (n == 1) { 
        fp = d[0][1] / h;
        fpp = 0;
    }

    fout << "\nf'(p) = " << fp << endl;
    fout << "f''(p) = " << fpp << endl;
}

int main() {
    ifstream fin("input.txt");
    ofstream fout("output.txt");

    double a, b, p;
    int n;

    fin >> a >> b >> n >> p;

    forwardDifferentiation(a, b, n, p, fout);

    fin.close();
    fout.close();
    return 0;
}
