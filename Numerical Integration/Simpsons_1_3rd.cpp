#include <bits/stdc++.h>
using namespace std;

#define vd vector<double>

ofstream fout;

double func(double x)
{
    return sqrt(x);
}

void simpsonOneThird(double a, double b, double h)
{
    fout<<"Simpson's 1/3:\n";
    int n = (b - a) / h;
    vd x(n + 1), fx(n + 1);

    for (int i = 0; i <= n; i++)
    {
        x[i] = a + i * h;
        fx[i] = func(x[i]);
    }

    double res = 0;
    for (int i = 0; i <= n; i++)
    {
        if (i == 0 || i == n)
            res += fx[i];
        else if (i % 2 != 0)
            res += 4 * fx[i];
        else
            res += 2 * fx[i];
    }

    fout <<"Result = "<<setprecision(3)<< res * (h / 3);
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

    simpsonOneThird(a, b, h);
    cout<<"Output in output.txt";
    fin.close();
    fout.close();
    return 0;
}
