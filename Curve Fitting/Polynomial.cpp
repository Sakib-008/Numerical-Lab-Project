#include <bits/stdc++.h>

using namespace std;

void printPolynomial(ofstream &out, vector<double> coeff, int m)
{
  out << "y = ";
  bool firstTerm = true;

  for (int i = 0; i < m; i++)
  {
    if (coeff[i] == 0)
      continue;

    if (!firstTerm)
      out << " + ";

    if (i == 0)
    {
      out << coeff[i];
    }
    else if (i == 1)
    {
      if (coeff[i] == 1)
        out << "x";
      else
        out << coeff[i] << "x";
    }
    else
    {
      if (coeff[i] == 1)
        out << "x^" << i;
      else
        out << coeff[i] << "x^" << i;
    }
    firstTerm = false;
  }
  out << endl;
}

void Polynomial(ofstream &out, vector<double> x, vector<double> y, int m, double newX, int n)
{
  int size = m + 1;
  vector<vector<double>> A(size, vector<double>(size, 0));
  vector<double> B(size, 0);

  for (int i = 0; i < size; i++)
  {
    for (int j = 0; j < size; j++)
    {
      for (int k = 0; k < n; k++)
      {
        A[i][j] += pow(x[k], i + j);
      }
    }
    for (int k = 0; k < n; k++)
    {
      B[i] += y[k] * pow(x[k], i);
    }
  }

  vector<double> coeff(size, 0);
  for (int i = 0; i < size; i++)
  {
    int maxi = i;
    for (int k = i + 1; k < size; k++)
    {
      if (fabs(A[k][i]) > fabs(A[maxi][i]))
      {
        maxi = k;
      }
    }
    if (A[maxi][i] == 0)
    {
      out << "Matrix is singular." << endl;
      return;
    }
    if (maxi != i)
    {
      swap(A[i], A[maxi]);
      swap(B[i], B[maxi]);
    }

    for (int j = i + 1; j < size; j++)
    {
      double factor = A[j][i] / A[i][i];
      for (int k = i; k < size; k++)
      {
        A[j][k] -= factor * A[i][k];
      }
      B[j] -= factor * B[i];
    }
  }
  for (int i = size - 1; i >= 0; i--)
  {
    coeff[i] = B[i];
    for (int j = i + 1; j < size; j++)
    {
      coeff[i] -= A[i][j] * coeff[j];
    }
    coeff[i] /= A[i][i];
  }

  printPolynomial(out, coeff, m);

  double newY = 0.0;
  for (int i = 0; i < size; i++)
  {
    newY += coeff[i] * pow(newX, i);
  }
  out << "For x = " << newX << ", predicted value of y is: " << newY << endl;
}

int main(void)
{
  ifstream in("PolynomialEquationInput.txt");
  ofstream out("PolynomialEquationOutput.txt");

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

  int m;
  in >> m;

  double newX;
  in >> newX;

  Polynomial(out, x, y, m, newX, n);

  in.close();
  out.close();
  return 0;
}