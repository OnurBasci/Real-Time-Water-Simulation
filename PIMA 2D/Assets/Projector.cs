using UnityEngine;
using System.Collections.Generic;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;
using System;

public class Projector
{
    int N = 16;
    public DenseMatrix u, v, p;
    double dt = 0.01;
    float density = 1f;
    const int SOLID = 0, FLUID = 1, EMPTY = 2;
    DenseMatrix types;

    public Projector()
    {
        // Initialize matrices
        u = DenseMatrix.Create(N + 1, N, 0.0);
        v = DenseMatrix.Create(N, N + 1, 0.0);
        p = DenseMatrix.Create(N, N, 0.0);

        // Initialize types matrix
        types = DenseMatrix.Create(N, N, EMPTY);
        for (int i = 0; i < N; i++)
        {
            types[0, i] = SOLID;
            types[i, 0] = SOLID;
            types[N - 1, i] = SOLID;
            types[i, N - 1] = SOLID;
        }

        RandomizeVectorField();
    }

    public Projector(DenseMatrix u, DenseMatrix v, DenseMatrix p, DenseMatrix t, int N, double dt, float density)
    {
        this.N = N;
        this.dt = dt;

        this.u = u;
        this.v = v;
        this.p = p;
        this.types = t;
        this.density = density;
    }

    private void RandomizeVectorField()
    {
        for (int i = 1; i < N - 1; i++)
        {
            for (int j = 1; j < N / 2; j++)
            {
                types[i, j] = FLUID;
                u[i, j] = UnityEngine.Random.value - 0.5;
                u[i + 1, j] = UnityEngine.Random.value - 0.5;
                v[i, j] = UnityEngine.Random.value - 0.5;
                v[i, j + 1] = UnityEngine.Random.value - 0.5;
            }
        }
    }

    public void Project()
    {
        //analyse fluid cells only
        List<Tuple<int, int>> fluidCells = new List<Tuple<int, int>>();
        Dictionary<Tuple<int, int>, int> fluidDict = new Dictionary<Tuple<int, int>, int>();

        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                if (types[i, j] == FLUID)
                {
                    fluidCells.Add(Tuple.Create(i, j));
                    fluidDict[Tuple.Create(i, j)] = fluidCells.Count - 1;
                }
            }
        }

        int nFluid = fluidCells.Count;
        var A = Matrix<double>.Build.Dense(nFluid, nFluid);
        var b = Vector<double>.Build.Dense(nFluid);

        for (int idx = 0; idx < nFluid; idx++)
        {
            int i = fluidCells[idx].Item1;
            int j = fluidCells[idx].Item2;

            if (types[i - 1, j] == FLUID)
                A[idx, fluidDict[Tuple.Create(i - 1, j)]] = -1;
            if (types[i, j - 1] == FLUID)
                A[idx, fluidDict[Tuple.Create(i, j - 1)]] = -1;
            if (types[i + 1, j] == FLUID)
                A[idx, fluidDict[Tuple.Create(i + 1, j)]] = -1;
            if (types[i, j + 1] == FLUID)
                A[idx, fluidDict[Tuple.Create(i, j + 1)]] = -1;

            int nonSolid = 0;
            double D = 0;

            //we suppose that solid's velocity is 0

            if (types[i - 1, j] != SOLID)
            {
                D -= u[i, j];
                nonSolid++;
            }
            if (types[i, j - 1] != SOLID)
            {
                D -= v[i, j];
                nonSolid++;
            }
            if (types[i + 1, j] != SOLID)
            {
                D += u[i + 1, j];
                nonSolid++;
            }
            if (types[i, j + 1] != SOLID)
            {
                D += v[i, j + 1];
                nonSolid++;
            }

            A[idx, idx] = nonSolid;
            b[idx] = -(this.density / dt) * D;
        }

        //Debug.Log(string.Join(", ", b.L2Norm()));

        // Print eigenvalues
        //var eigenvalues = A.Evd().EigenValues;
        //Console.WriteLine("eig " + eigenvalues);

        // Solve linear system
        var x = A.Solve(b);

        // Convert x to 2D array
        double[,] p = new double[N, N];
        foreach (var cell in fluidCells)
        {
            int i = cell.Item1;
            int j = cell.Item2;
            p[i, j] = x[fluidDict[Tuple.Create(i, j)]];
        }

        // Update u2
        double[,] u2 = new double[N, N];
        for (int i = 0; i < N - 1; i++)
        {
            for (int j = 0; j < N; j++)
            {
                if (types[i, j] != FLUID && types[i + 1, j] != FLUID)
                    continue;

                double rightP = 0, leftP = 0;

                if (types[i + 1, j] == FLUID)
                    rightP = x[fluidDict[Tuple.Create(i + 1, j)]];
                if (types[i + 1, j] == EMPTY)
                    rightP = 0;
                if (types[i + 1, j] == SOLID)
                {
                    u2[i + 1, j] = 0;
                    continue;
                }

                if (types[i, j] == FLUID)
                    leftP = x[fluidDict[Tuple.Create(i, j)]];
                if (types[i, j] == EMPTY)
                    leftP = 0;
                if (types[i, j] == SOLID)
                {
                    u2[i + 1, j] = 0;
                    continue;
                }

                u2[i + 1, j] += dt * (rightP - leftP);
            }
        }

        // Update v2
        double[,] v2 = new double[N, N];
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N - 1; j++)
            {
                if (types[i, j] != FLUID && types[i, j + 1] != FLUID)
                    continue;

                double topP = 0, bottomP = 0;

                if (types[i, j + 1] == FLUID)
                    topP = x[fluidDict[Tuple.Create(i, j + 1)]];
                if (types[i, j + 1] == EMPTY)
                    topP = 0;
                if (types[i, j + 1] == SOLID)
                {
                    v2[i, j + 1] = 0;
                    continue;
                }

                if (types[i, j] == FLUID)
                    bottomP = x[fluidDict[Tuple.Create(i, j)]];
                if (types[i, j] == EMPTY)
                    bottomP = 0;
                if (types[i, j] == SOLID)
                {
                    v2[i, j + 1] = 0;
                    continue;
                }

                v2[i, j + 1] += dt * (topP - bottomP);
            }
        }


        // Update u and v
        u = DenseMatrix.OfArray(u2);
        v = DenseMatrix.OfArray(v2);
        this.p = DenseMatrix.OfArray(p);
    }
}
