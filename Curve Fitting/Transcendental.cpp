#include <bits/stdc++.h>

using namespace std;

// The transcendental equation format is: y=ax^b
// lny = lna + blnx (Taking ln on both sides)
// Y = A + bX
// Here, Y = lny, X = lnx, A = lna
void Transcendental(ofstream &out, vector<double> x, vector<double> y, int newX)
{
  int n = x.size();
  double Sx = 0, Sy = 0, Sxx = 0, Sxy = 0;
  for (int i = 0; i < n; i++)
  {
    Sx += log(x[i]);
    Sy += log(y[i]);
    Sxx += log(x[i]) * log(x[i]);
    Sxy += log(x[i]) * log(y[i]);
  }
  double b = (n * Sxy - Sx * Sy) / (n * Sxx - Sx * Sx);

  double A = (Sy - b * Sx) / n;
  double a = exp(A);

  out << "The Transcendental Equation is: " << endl;
  out << "y = " << a << "x^" << b << endl;

  double newY = a * pow(newX, b);
  out << "For x = " << newX << ", the predicted value of y is: ";
  out << newY << endl;
}

int main(void)
{
  ifstream in("TranscendentalEquationInput.txt");
  ofstream out("TranscendentalEquationOutput.txt");

  int n;
  in >> n;
  vector<double> x(n), y(n);

  for (int i = 0; i < n; i++)
  {
    in >> x[i] >> y[i];
  }

  out << setprecision(2) << fixed;

  out << "Data points entered:" << endl;
  for (int i = 0; i < n; i++)
  {
    out << "(" << x[i] << ", " << y[i] << ")\n";
  }

  int newX;
  in >> newX;
  Transcendental(out, x, y, newX);

  in.close();
  out.close();
  return 0;
}