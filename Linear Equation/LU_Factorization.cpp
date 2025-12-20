// LU Decomposition implemented by 2207008

#include <bits/stdc++.h>
using namespace std;

int main()
{
    ifstream in("LU_input.txt");
    ofstream out("LU_output.txt");

    if (!in || !out)
    {
        cout << "Error opening input/output file!" << endl;
        return 1;
    }

    out << fixed << setprecision(3);
    out << "LU Factorization method : " << endl;

    int t; // Number of test cases
    in >> t;
    for (int cs = 1; cs <= t; cs++)
    {
        int n;
        in >> n; // Number of equations.
        vector<vector<double>> A(n, vector<double>(n)), L(n, vector<double>(n, 0)), U(n, vector<double>(n, 0));
        vector<double> B(n), X(n), Y(n);

        // Entering the augmented matrix.
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                in >> A[i][j];
            }
            in >> B[i];
        }
        for (int i = 0; i < n; i++)
        {
            for (int j = i; j < n; j++)
            {
                double sum = 0;
                for (int k = 0; k < i; k++)
                {
                    sum += L[i][k] * U[k][j];
                }
                U[i][j] = A[i][j] - sum;
            }
            for (int j = i; j < n; j++)
            {
                if (i == j)
                    L[i][j] = 1;
                else
                {
                    double sum = 0;
                    for (int k = 0; k < i; k++)
                    {
                        sum += L[j][k] * U[k][i];
                    }
                    L[j][i] = (A[j][i] - sum) / U[i][i];
                }
            }
        }
        out << "Output for case " << cs << " : " << endl;
        out << "\nLower Triangular Matrix (L) : " << endl;
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                out << L[i][j] << " ";
            }
            out << endl;
        }
        out << "\nUpper Triangular Matrix (U) : " << endl;
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                out << U[i][j] << " ";
            }
            out << endl;
        }
        for (int i = 0; i < n; i++)
        {
            double sum = 0;
            for (int j = 0; j < i; j++)
            {
                sum += L[i][j] * Y[j];
            }
            Y[i] = B[i] - sum;
        }
        bool soln = true;
        for (int i = 0; i < n; i++)
        {
            if (fabs(U[i][i]) < 1e-9 && fabs(Y[i]) > 1e-9)
            {
                out << "\nNo Solution!" << endl;
                out << "\n------------------------------\n"
                    << endl;
                soln = false;
                break;
            }
            else if (fabs(U[i][i]) < 1e-9 && fabs(Y[i]) < 1e-9)
            {
                out << "\nInfinite Solution!" << endl;
                out << "\n------------------------------\n"
                    << endl;
                soln = false;
                break;
            }
        }
        if (soln == false)
            continue;
        for (int i = n - 1; i >= 0; i--)
        {
            double sum = 0;
            for (int j = i + 1; j < n; j++)
            {
                sum += U[i][j] * X[j];
            }
            X[i] = (Y[i] - sum) / U[i][i];
        }
        out << "\nSystem has unique solution!" << endl;
        out << "\nSolution : " << endl;
        for (int i = 0; i < n; i++)
        {
            out << "x" << i + 1 << " = " << X[i] << endl;
        }
        out << "\n------------------------------\n"
            << endl;
    }
    in.close();
    out.close();
    return 0;
}
