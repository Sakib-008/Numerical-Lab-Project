// Gauss Elimination Method implemented by 2207009
#include <bits/stdc++.h>
using namespace std;

void printMatrix(ofstream &out, vector<vector<double>> mat)
{
  int n = mat.size();
  int m = mat[0].size();
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < m; j++)
      out << mat[i][j] << " ";
    out << endl;
  }
}

void gaussElimination(ofstream &out, vector<vector<double>> A, int n)
{
  int flag = 0;

  for (int i = 0; i < n; i++)
  {
    if (A[i][i] == 0)
    {
      for (int r = i + 1; r < n; r++)
      {
        if (A[r][i] != 0)
        {
          swap(A[i], A[r]);
          break;
        }
      }
      if (A[i][i] == 0)
      {
        flag = 1;
        break;
      }
    }

    for (int j = i + 1; j < n; j++)
    {
      double ratio = A[j][i] / A[i][i];
      for (int k = 0; k < n + 1; k++)
      {
        A[j][k] = A[j][k] - ratio * A[i][k];
      }
    }
  }

  if (flag == 1)
  {
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

      if (allZero && A[i][n] != 0)
        noSolution = true;
      if (allZero && A[i][n] == 0)
        infiniteSolution = true;
    }
    if (noSolution)
    {
      out << "No solution" << endl;
    }
    else if (infiniteSolution)
    {
      out << "Infinite solutions" << endl;
    }
  }
  else
  {
    out << "The Row Echelon Form is: " << endl;
    printMatrix(out, A);
    vector<double> X(n);

    for (int i = n - 1; i >= 0; i--)
    {
      double sum = A[i][n];
      for (int j = i + 1; j < n; j++)
      {
        sum -= A[i][j] * X[j];
      }
      X[i] = sum / A[i][i];
    }

    out << "The unique solution is: " << endl;
    for (int i = 0; i < n; i++)
      out << "x" << i + 1 << " = " << X[i] << (i == n - 1 ? "" : " , ");
    out << endl;
  }
}

int main(void)
{
  ifstream in("GaussEliminationInput.txt");
  ofstream out("GaussEliminationOutput.txt");

  int T;
  in >> T;

  while (T--)
  {
    int n;
    in >> n;

    vector<vector<double>> A(n, vector<double>(n + 1));

    for (int i = 0; i < n; i++)
    {
      for (int j = 0; j < n; j++)
        in >> A[i][j];
      in >> A[i][n];
    }

    out << setprecision(3) << fixed;

    out << "The system of equations are: " << endl;
    for (int i = 0; i < n; i++)
    {
      for (int j = 0; j < n; j++)
      {
        if (j > 0)
          out << (A[i][j] < 0 ? " - " : " + ");

        out << abs(A[i][j]) << "x" << j + 1;
      }
      out << " = " << A[i][n] << endl;
    }
    gaussElimination(out, A, n);
  }

  in.close();
  out.close();

  return 0;
}