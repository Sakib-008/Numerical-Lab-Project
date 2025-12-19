#include <bits/stdc++.h>
using namespace std;

#define vd vector<double>

ofstream fout;

double func(double x)
{
    return sqrt(x);
}


void simpsonThreeEight(double a, double b, double h)
{
    fout<<"Simpson's 3/8:\n";
    int n = round((b - a) / h);
    if (n % 3 != 0)
    {
        fout << "Error: n must be multiple of 3 for Simpson 3/8 rule\n";
        return;
    }

    h = (b - a) / n;

    vd x(n + 1), fx(n + 1);
    for (int i = 0; i <= n; i++)
    {
        x[i] = a + i * h;
        fx[i] = func(x[i]);
    }

    double res = fx[0] + fx[n];
    for (int i = 1; i < n; i++)
    {
        if (i % 3 == 0)
            res += 2 * fx[i];
        else
            res += 3 * fx[i];
    }

    fout <<setprecision(3)<< res * (3 * h / 8);
}

int main()
{
    // input format
    //a b h
    ifstream fin("input.txt");
    fout.open("output.txt");

    double a, b, n;
    fin >> a >> b >> n;
    double h = (b-a)/n;

    simpsonThreeEight(a, b, h);
    cout<<"Output in output.txt";
    fin.close();
    fout.close();
    return 0;
}
