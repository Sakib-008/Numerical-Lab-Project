// Matrix Inversion Method implemented by 2207009
#include <bits/stdc++.h>
using namespace std;

vector<vector<double>> cofactorMatrix(vector<vector<double>> A, int p, int q)
{
  int n = A.size();
  vector<vector<double>> temp;
  temp.reserve(n - 1);
  for (int i = 0; i < n; i++)
  {
    if (i == p)
      continue;
    vector<double> row;
    row.reserve(n - 1);
    for (int j = 0; j < n; j++)
    {
      if (j == q)
        continue;
      row.push_back(A[i][j]);
    }
    temp.push_back(row);
  }
  return temp;
}

double determinantMatrix(vector<vector<double>> A)
{
  int n = A.size();
  if (n == 1)
    return A[0][0];
  if (n == 2)
    return A[0][0] * A[1][1] - A[0][1] * A[1][0];

  double detA = 0.0;

  for (int j = 0; j < n; j++)
  {
    vector<vector<double>> cofactor = cofactorMatrix(A, 0, j);
    double value = determinantMatrix(cofactor);
    if (j % 2 == 1)
      value = -value;
    detA += value * A[0][j];
  }
  return detA;
}

vector<vector<double>> transposeMatrix(vector<vector<double>> A)
{
  int n = A.size();
  vector<vector<double>> T(n, vector<double>(n));
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < n; j++)
      T[j][i] = A[i][j];
  }
  return T;
}

vector<vector<double>> adjointMatrix(vector<vector<double>> A)
{
  int n = A.size();
  if (n == 1)
  {
    A[0][0] = 1;
    return A;
  };
  vector<vector<double>> adjoint(n, vector<double>(n));
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < n; j++)
    {
      vector<vector<double>> cofactor = cofactorMatrix(A, i, j);
      double value = determinantMatrix(cofactor);
      if ((i + j) % 2 == 1)
        value = -value;
      adjoint[i][j] = value;
    }
  }
  return transposeMatrix(adjoint);
}

vector<vector<double>> inverseMatrix(vector<vector<double>> A)
{
  double detA = determinantMatrix(A);
  if (detA == 0)
    throw runtime_error("Singular Matrix");
  int n = A.size();
  vector<vector<double>> adjoint = adjointMatrix(A);
  vector<vector<double>> inverse(n, vector<double>(n));
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < n; j++)
      inverse[i][j] = adjoint[i][j] / detA;
  }
  return inverse;
}

vector<double> calculateSolution(vector<vector<double>> invA, vector<double> B)
{
  int n = B.size();
  vector<double> X(n, 0);
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < n; j++)
      X[i] += invA[i][j] * B[j];
  }
  return X;
}

void printMat(ofstream &out, vector<vector<double>> A)
{
  int n = A.size();
  int m = A[0].size();
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < n; j++)
      out << A[i][j] << " ";
    out << endl;
  }
}

int main(void)
{
  ifstream in("MatrixInversionInput.txt");
  ofstream out("MatrixInversionOutput.txt");

  int T;
  in >> T;

  while (T--)
  {

    int n;
    in >> n;

    vector<vector<double>> A(n, vector<double>(n));
    vector<double> B(n);

    for (int i = 0; i < n; i++)
    {
      for (int j = 0; j < n; j++)
        in >> A[i][j];
      in >> B[i];
    }

    out << setprecision(3) << fixed;

    out << "The system of equations are: " << endl;
    for (int i = 0; i < n; i++)
    {
      for (int j = 0; j < n; j++)
      {
        if (j > 0)
        {
          out << (A[i][j] < 0 ? " - " : " + ");
        }
        out << abs(A[i][j]) << "x" << j + 1;
      }
      out << " = " << B[i] << endl;
    }

    double detA = determinantMatrix(A);

    if (detA > 0)
    {
      vector<vector<double>> transposeMat = adjointMatrix(A);
      out << "The transpose matrix is: " << endl;
      printMat(out, transposeMat);
      vector<vector<double>> invA = inverseMatrix(A);
      out << "The final inverse matrix is: " << endl;
      printMat(out, invA);
      vector<double> X = calculateSolution(invA, B);
      out << "The unique solution is: " << endl;
      for (int i = 0; i < n; i++)
        out << "x" << i + 1 << " = " << X[i] << (i == n - 1 ? "" : " , ");
      out << endl;
    }
    else
    {
      for (int i = 0; i < n; i++)
      {
        if (A[i][i] == 0)
        {
          for (int r = i + 1; r < n; r++)
          {
            if (A[r][i] != 0)
            {
              swap(A[i], A[r]);
              swap(B[i], B[r]);
              break;
            }
          }
          if (A[i][i] == 0)
            break;
        }

        for (int j = i + 1; j < n; j++)
        {
          double ratio = A[j][i] / A[i][i];
          for (int k = 0; k < n; k++)
            A[j][k] = A[j][k] - ratio * A[i][k];
          B[j] = B[j] - ratio * B[i];
        }
      }

      bool noSolution = false;
      bool infiniteSolution = false;
      for (int i = 0; i < n; i++)
      {
        bool allZero = true;
        for (int j = 0; j < n; j++)
        {
          if (A[i][j] > 0)
          {
            allZero = false;
            break;
          }
        }

        if (allZero && B[i] != 0)
          noSolution = true;
        if (allZero && B[i] == 0)
          infiniteSolution = true;
      }

      if (noSolution)
        out << "No solution" << endl;
      else if (infiniteSolution)
        out << "Infinite solutions" << endl;
    }
  }
  in.close();
  out.close();
  return 0;
}