using System;
using System.Collections.Generic;
using System.Linq;
using Unity.Mathematics;
using UnityEngine;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;

public class fluidShowerV2 : MonoBehaviour
{
    public GameObject squarePrefab; // Prefab for the grid square
    public float squareSize = 1f; // Size of each square

    public float gravity = -9.81f;
    public float density = 1000f;
    public int numX = 100;
    public int numY = 100;
    public float h = 1;
    public FieldType currentField;
    public bool activateHorizontalFore = false;
    public bool addCircleToMiddle = false;
    public float horizontalForce = 40;
    public bool printVolume = false;
    public bool printGrid = false;

    private float deltaTime;

    private double[,] bufferZeroArray;

    public enum FieldType
    {
        MFIELD,
        PFIELD,
        UFIELD,
        VFIELD
    }

    [NonSerialized] public Dictionary<GameObject, int2> blockInfos;

    public const int UFIELD = 0;
    public const int VFIELD = 1;
    public const int SFIELD = 2; //smoke field

    private GameObject[,] blockObjectsArray;

    [NonSerialized] public Fluid fluid;

    public void Start()
    {
        blockObjectsArray = new GameObject[numX, numY];

        fluid = new Fluid(density, numX, numY, h, this);
        GenerateGrid(fluid);
        //fluid.setUp();

        if(addCircleToMiddle)
        {
            fluid.addCircleFluidToMiddle();
        }


        bufferZeroArray = new double[numX * numY, numX * numY];

        // Initialize each element of the array to 0
        for (int i = 0; i < numX * numY; i++)
        {
            for (int j = 0; j < numX * numY; j++)
            {
                bufferZeroArray[i, j] = 0.0;
            }
        }
    }

    public void Update()
    {
        //get the delta time
        deltaTime += (Time.deltaTime - deltaTime) * 0.1f;

        fluid.Simulate(deltaTime, gravity);
        simulateGrid(fluid);

        if(printVolume)
        {
            fluid.printVolume();
        }

        if(printGrid)
        {
            fluid.printGrid();
        }
    }


    void GenerateGrid(Fluid fluid)
    {
        blockInfos = new Dictionary<GameObject, int2>();

        //if (max_pressure <= 0) max_pressure = gridArray.Min();
        for (int i = 0; i < numY; i++)
        {
            for (int j = 0; j < numX; j++)
            {
                Vector3 spawnPosition = new Vector3(j * squareSize, i * squareSize, 0);
                GameObject square = Instantiate(squarePrefab, spawnPosition, Quaternion.identity);
                square.transform.localScale = new Vector3(squareSize, squareSize, 1f);

                blockInfos.Add(square, new int2(j, i));

                blockObjectsArray[j, i] = square;
            }
        }
    }

    public void simulateGrid(Fluid fluid)
    {
        double[] gridArray = fluid.m;
        switch (currentField)
        {
            case FieldType.MFIELD:
                gridArray = fluid.m;
                break;
            case FieldType.PFIELD:
                gridArray = fluid.p;
                break;
            case FieldType.UFIELD:
                gridArray = fluid.u;
                break;
            case FieldType.VFIELD:
                gridArray = fluid.v;
                break;
        }

        //get maximum pressure
        double max_pressure = gridArray.Max();
        double min_pressure = gridArray.Min();


        for (int i = 0; i < numY; i++)
        {
            for (int j = 0; j < numX; j++)
            {

                // Set color based on gridArray value
                double c = (gridArray[j * numY + i] - min_pressure) / (max_pressure - min_pressure);

                UnityEngine.Color color = new UnityEngine.Color((float)c, (float)c, (float)c);

                blockObjectsArray[j, i].GetComponent<SpriteRenderer>().color = color;
            }
        }
    }


    public class Fluid
    {
        public fluidShowerV2 fluidShower;
        public float density;
        public int numX;
        public int numY;
        public int numCells;
        public float h;
        public double[] u;
        public double[] v;
        public double[] newU;
        public double[] newV;
        public double[] p;
        public int[] s;
        public double[] m;  // m for density (mase)
        public double[] newM;

        public Fluid(float density, int numX, int numY, float h, fluidShowerV2 fluidShower)
        {
            this.density = density;
            this.numX = numX;
            this.numY = numY;
            this.numCells = this.numX * this.numY;
            this.h = h;
            this.u = new double[this.numCells];
            this.v = new double[this.numCells];
            this.newU = new double[this.numCells];
            this.newV = new double[this.numCells];
            this.p = new double[this.numCells];
            this.s = new int[this.numCells];
            this.m = new double[this.numCells];
            this.newM = new double[this.numCells];
            this.fluidShower = fluidShower;

            // Initialize m and s array with 1.0
            for (int i = 0; i < numX; i++)
            {
                for (int j = 0; j < numY; j++)
                {
                    //this.m[i] = 1.0f;
                    this.s[i * numY + j] = 1;
                    if (i == 0 || i == numX - 1 || j == 0 || j == numY - 1)
                        s[i * numY + j] = 0;	// make the border solid
                }
            }

            //m[5 * numY + 10] = 1;

        }
        private void Integrate(float dt, float gravity)
        {
            int n = this.numY;
            for (int i = 1; i < this.numX; i++)
            {
                for (int j = 1; j < this.numY - 1; j++)
                {
                    if (this.s[i * n + j] != 0.0f && this.s[i * n + j - 1] != 0.0f)
                    {
                        this.v[i * n + j] += gravity * dt;
                    }
                }
            }
        }

        private void SolveIncompressibility(float dt)
        {
            DenseMatrix u = DenseMatrix.Create(numX, numY, 0);
            DenseMatrix v = DenseMatrix.Create(numX, numY, 0);
            DenseMatrix p = DenseMatrix.Create(numX, numY, 0);

            for(int i = 0; i < numX; i ++)
            {
                for(int j = 0; j < numY; j++)
                {
                    u[i,j] = this.u[i*numY + j];
                    v[i,j] = this.v[i*numY + j];
                    p[i,j] = this.p[i*numY + j];
                }
            }

            Projector projector = new Projector(u, v, p, numX, dt, density);
            projector.Project();

            for (int i = 0; i < numX; i++)
            {
                for (int j = 0; j < numY; j++)
                {
                    this.u[i * numY + j] = projector.u[i,j];
                    this.v[i * numY + j] = projector.v[i,j];
                    this.p[i * numY + j] = projector.p[i,j];
                }
            }
        }

        private void updatePressureField(float dt)
        {
            //In this method we solve the linear equation Ap = d for p where A is the pressure coefficients, p is the pressure
            //values and d is the divergence of a block

            //calculate A and d
            //var A = new DenseMatrix(numX*numY, numX*numY);

            var A = Matrix<double>.Build.DenseOfArray(fluidShower.bufferZeroArray);
            Vector<double> divergences = Vector<double>.Build.Dense(numX*numY);

            //we fill the diagonal coefficient with the number of neighboors
            var n = numY;
            for (int i = 1; i < numX - 1;  i ++)
            {
                for(int j = 1; j < numY - 1; j ++)
                {
                    double sx0 = this.s[(i - 1) * n + j];
                    double sx1 = this.s[(i + 1) * n + j];
                    double sy0 = this.s[i * n + j - 1];
                    double sy1 = this.s[i * n + j + 1];
                    double nbNeighboors = sx0 + sx1 + sy0 + sy1;

                    //we fill the diagonal coefficient with the number of neighboors
                    A[i * numY + j, i * numY + j] = nbNeighboors;

                    //others are filled with -1: pi+1,j, pi,j+1, pi-1,j , pi,j-1
                    A[i * numY + j, (i + 1) * numY + j] = -1;
                    A[i * numY + j, i * numY + j + 1] = -1;
                    A[i * numY + j, (i-1) * numY + j] = -1;
                    A[i * numY + j, i * numY + j - 1] = -1;

                    double div = (u[(i + 1) * n + j] - u[i * n + j] + v[i * n + j + 1] - v[i * n + j])/h;

                    divergences[i * numY + j] = div;
                }
            }

            double divCosntant = density * math.pow(h, 2) / dt;
            divergences *= divCosntant;

            Debug.Log("Before solving");

            Vector<double> pressureVector = A.Solve(divergences);

            Debug.Log("After solving");

            //update the pressure field
            //p = pressureVector.ToColumnMatrix().ToArray2D();

            for (int i = 0; i < numX; i++)
            {
                for(int j = 0; j < numY; j++)
                {
                    p[i* numY + j] = (float) pressureVector[i * numY + j]; 
                }
            }
        }

        private void updatePressureFieldV2(float dt)
        {
            List<Tuple<int, int>> fluidCells = new List<Tuple<int, int>>();
            Dictionary<Tuple<int, int>, int> fluidDict = new Dictionary<Tuple<int, int>, int>();

            for (int i = 0; i < numX; i++)
            {
                for (int j = 0; j < numY; j++)
                { 
                    if (s[i * numY + j] == 1)
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
                //Calculate the pressure coefficients
                int i = fluidCells[idx].Item1;
                int j = fluidCells[idx].Item2;

                if (i == 0 || i == numX - 1 || j == 0 || j == numY - 1)
                    continue;

                if (s[(i - 1)*numY + j] == 1)
                    A[idx, fluidDict[Tuple.Create(i - 1, j)]] = -1;
                if (s[i * numY + j - 1] == 1)
                    A[idx, fluidDict[Tuple.Create(i, j - 1)]] = -1;
                if (s[(i + 1) * numY + j] == 1)
                    A[idx, fluidDict[Tuple.Create(i + 1, j)]] = -1;
                if (s[i * numY + j + 1] == 1)
                    A[idx, fluidDict[Tuple.Create(i, j + 1)]] = -1;

                //calculate the diagonal coefficient
                int nonSolid = 0;
                double D = 0;

                if (s[(i - 1) * numY + j] != 0)
                {
                    D -= u[i * numY + j];
                    nonSolid++;
                }
                if (s[i * numY + j - 1] != 0)
                {
                    D -= v[i * numY + j];
                    nonSolid++;
                }
                if (s[(i + 1) * numY + j] != 0)
                {
                    D += u[(i + 1) * numY + j];
                    nonSolid++;
                }
                if (s[i * numY + j + 1] != 0)
                {
                    D += v[i * numY + j + 1];
                    nonSolid++;
                }

                A[idx, idx] = nonSolid;
                b[idx] = -(1 / dt) * D;
            }

            // Solve linear system
            var x = A.Solve(b);

            Debug.Log(string.Join(", ", x));

            // Convert x to 2D array
            foreach (var cell in fluidCells)
            {
                int i = cell.Item1;
                int j = cell.Item2;
                p[i*numY+j] = x[fluidDict[Tuple.Create(i, j)]];
            }

            // Update u2
            double[] u2 = new double[numCells];
            for (int i = 0; i < numX - 1; i++)
            {
                for (int j = 0; j < numY; j++)
                {
                    if (s[i*numX+ j] != 1 && s[(i + 1) * numY + j] != 1)
                        continue;

                    double rightP = 0, leftP = 0;

                    if (s[(i + 1) * numY + j] == 1)
                        rightP = x[fluidDict[Tuple.Create(i + 1, j)]];
                    if (s[(i + 1) * numY + j] == 2) //If it is empty to add
                        rightP = 0;
                    if (s[(i + 1) * numY + j] == 0)
                    {
                        u2[(i + 1) * numY + j] = 0;
                        continue;
                    }

                    if (s[(i * numY) + j] == 1)
                        leftP = x[fluidDict[Tuple.Create(i, j)]];
                    if (s[(i * numY) + j] == 2) //TO ADD
                        leftP = 0;
                    if (s[i * numY + j] == 0)
                    {
                        u2[(i + 1) * numY + j] = 0;
                        continue;
                    }

                    u2[(i + 1) * numY + j] -= dt * (rightP - leftP);
                }
            }

            // Update v2
            double[] v2 = new double[numCells];
            for (int i = 0; i < numX; i++)
            {
                for (int j = 0; j < numY - 1; j++)
                {
                    if (s[i* numY + j] != 1 && s[i*numY +j + 1] != 1)
                        continue;

                    double topP = 0, bottomP = 0;

                    if (s[i * numY + j + 1] == 1)
                        topP = x[fluidDict[Tuple.Create(i, j + 1)]];
                    if (s[i * numY + j + 1] == 2) //TO ADD
                        topP = 0;
                    if (s[i * numY + j + 1] == 0)
                    {
                        v2[i * numY + j + 1] = 0;
                        continue;
                    }

                    if (s[i * numY + j] == 1)
                        bottomP = x[fluidDict[Tuple.Create(i, j)]];
                    if (s[i * numY + j] == 2)
                        bottomP = 0;
                    if (s[i * numY + j] == 0)
                    {
                        v2[i * numY + j + 1] = 0;
                        continue;
                    }

                    v2[i * numY + j + 1] -= dt * (topP - bottomP);
                }
            }

            // Update u and v
            u = u2;
            v = v2;
        }

        private void Extrapolate()
        {
            int n = this.numY;

            // Extrapolate u values at top and bottom boundaries
            for (int i = 0; i < this.numX; i++)
            {
                this.u[i * n + 0] = this.u[i * n + 1];
                this.u[i * n + this.numY - 1] = this.u[i * n + this.numY - 2];
            }

            // Extrapolate v values at left and right boundaries
            for (int j = 0; j < this.numY; j++)
            {
                this.v[0 * n + j] = this.v[1 * n + j];
                this.v[(this.numX - 1) * n + j] = this.v[(this.numX - 2) * n + j];
            }
        }

        private void Extrapolate2()
        {
            int N = numX;

            // Define maximum integer value for marking unreachable cells
            int maxint = N * N;

            // Array to store distances to nearest FLUID cells
            int[,] d = new int[N + 1, N];

            // List to store cells to be processed
            List<Tuple<int, int>> W = new List<Tuple<int, int>>();

            // First pass for u
            for (int i = 0; i < N - 1; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    if (s[i*numY + j] == 1 || s[i*numY + j] == 1)
                    {
                        d[i + 1, j] = 0;
                    }
                }
            }

            // Initialize W list with cells having a FLUID neighbor
            for (int i = 0; i < N + 1; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    if (d[i, j] == 0)
                    {
                        continue;
                    }

                    foreach (var (i2, j2) in new[] { (i - 1, j), (i, j - 1), (i + 1, j), (i, j + 1) })
                    {
                        if (i2 >= 0 && i2 < N + 1 && j2 >= 0 && j2 < N)
                        {
                            if (d[i2, j2] == 0)
                            {
                                d[i, j] = 1;
                                W.Add(Tuple.Create(i, j));
                                break;
                            }
                        }
                    }
                }
            }

            // Perform the main extrapolation loop
            int t = 0;
            while (t < W.Count)
            {
                var (i, j) = W[t];
                double newU = 0;
                int count = 0;
                foreach (var (i2, j2) in new[] { (i - 1, j), (i, j - 1), (i + 1, j), (i, j + 1) })
                {
                    if (i2 >= 0 && i2 < N + 1 && j2 >= 0 && j2 < N)
                    {
                        //check if the cell is closer to boundary
                        if (d[i2, j2] < d[i, j])
                        {
                            newU += u[i2*numY + j2];
                            count++;
                        }
                    }
                }
                u[i*numY + j] = newU / count;

                foreach (var (i2, j2) in new[] { (i - 1, j), (i, j - 1), (i + 1, j), (i, j + 1) })
                {
                    if (i2 >= 0 && i2 < N + 1 && j2 >= 0 && j2 < N)
                    {
                        if (d[i2, j2] == maxint)
                        {
                            d[i2, j2] = d[i, j] + 1;
                            W.Add(Tuple.Create(i2, j2));
                        }
                    }
                }
                t++;
            }

            // Clear distance array for the next pass
            d = new int[N, N + 1];
            W.Clear();

            // Second pass for v
            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < N - 1; j++)
                {
                    if (s[i*numY+j] == 1 || s[i*numY+ j + 1] == 1)
                    {
                        d[i, j + 1] = 0;
                    }
                }
            }

            // Initialize W list with cells having a FLUID neighbor
            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < N + 1; j++)
                {
                    if (d[i, j] == 0)
                    {
                        continue;
                    }

                    foreach (var (i2, j2) in new[] { (i - 1, j), (i, j - 1), (i + 1, j), (i, j + 1) })
                    {
                        if (i2 >= 0 && i2 < N && j2 >= 0 && j2 < N + 1)
                        {
                            if (d[i2, j2] == 0)
                            {
                                d[i, j] = 1;
                                W.Add(Tuple.Create(i, j));
                                break;
                            }
                        }
                    }
                }
            }

            // Perform the main extrapolation loop for v
            t = 0;
            while (t < W.Count)
            {
                var (i, j) = W[t];
                double newV = 0;
                int count = 0;
                foreach (var (i2, j2) in new[] { (i - 1, j), (i, j - 1), (i + 1, j), (i, j + 1) })
                {
                    if (i2 >= 0 && i2 < N && j2 >= 0 && j2 < N + 1)
                    {
                        if (d[i2, j2] < d[i, j])
                        {
                            newV += v[i2*numY+ j2];
                            count++;
                        }
                    }
                }
                v[i*numY+j] = newV / count;

                foreach (var (i2, j2) in new[] { (i - 1, j), (i, j - 1), (i + 1, j), (i, j + 1) })
                {
                    if (i2 >= 0 && i2 < N && j2 >= 0 && j2 < N + 1)
                    {
                        if (d[i2, j2] == maxint)
                        {
                            d[i2, j2] = d[i, j] + 1;
                            W.Add(Tuple.Create(i2, j2));
                        }
                    }
                }
                t++;
            }
        }

        public double SampleField(double x, double y, int field)
        {
            //this function applies a biliniar interpolation to calculate the value of a point in a specified field
            int n = this.numY;
            float h = this.h;
            float h1 = 1.0f / h;
            float h2 = 0.5f * h;

            x = Mathf.Max(Mathf.Min((float)x, this.numX * h), h);
            y = Mathf.Max(Mathf.Min((float)y, this.numY * h), h);

            float dx = 0.0f;
            float dy = 0.0f;

            double[] f = null;

            switch (field)
            {
                case UFIELD:
                    f = this.u;
                    dy = h2;
                    break;
                case VFIELD:
                    f = this.v;
                    dx = h2;
                    break;
                case SFIELD:
                    f = this.m;
                    dx = h2;
                    dy = h2;
                    break;
            }

            int x0 = Mathf.Min(Mathf.FloorToInt((float)(x - dx) * h1), this.numX - 1);
            float tx = ((float)(x - dx) - x0 * h) * h1;
            int x1 = Mathf.Min(x0 + 1, this.numX - 1);

            int y0 = Mathf.Min(Mathf.FloorToInt((float)(y - dy) * h1), this.numY - 1);
            float ty = ((float)(y - dy) - y0 * h) * h1;
            int y1 = Mathf.Min(y0 + 1, this.numY - 1);

            float sx = 1.0f - tx;
            float sy = 1.0f - ty;

            double val = sx * sy * f[x0 * n + y0] +
                        tx * sy * f[x1 * n + y0] +
                        tx * ty * f[x1 * n + y1] +
                        sx * ty * f[x0 * n + y1];

            return val;
        }

        public double AvgU(int i, int j)
        {
            int n = this.numY;
            double u = (this.u[i * n + j - 1] + this.u[i * n + j] +
                       this.u[(i + 1) * n + j - 1] + this.u[(i + 1) * n + j]) * 0.25f;
            return u;
        }

        public double AvgV(int i, int j)
        {
            int n = this.numY;
            double v = (this.v[(i - 1) * n + j] + this.v[i * n + j] +
                       this.v[(i - 1) * n + j + 1] + this.v[i * n + j + 1]) * 0.25f;
            return v;
        }

        public void AdvectVel(float dt)
        {
            // Copy current velocities to temporary arrays
            this.newU = (double[])this.u.Clone();
            this.newV = (double[])this.v.Clone();

            int n = this.numY;
            float h = this.h;
            float h2 = 0.5f * h;

            for (int i = 1; i < this.numX; i++)
            {
                for (int j = 1; j < this.numY; j++)
                {
                    // u component
                    if (this.s[i * n + j] != 0.0f && this.s[(i - 1) * n + j] != 0.0f && j < this.numY - 1)
                    {
                        double x = i * h;
                        double y = j * h + h2;
                        double u = this.u[i * n + j];
                        double v = this.AvgV(i, j);
                        x = x - dt * u;
                        y = y - dt * v;
                        u = this.SampleField(x, y, UFIELD);
                        this.newU[i * n + j] = u;
                    }
                    // v component
                    if (this.s[i * n + j] != 0.0f && this.s[i * n + j - 1] != 0.0f && i < this.numX - 1)
                    {
                        double x = i * h + h2;
                        double y = j * h;
                        double u = this.AvgU(i, j);
                        double v = this.v[i * n + j];
                        x = x - dt * u;
                        y = y - dt * v;
                        v = this.SampleField(x, y, VFIELD);
                        this.newV[i * n + j] = v;
                    }
                }
            }

            // Update velocities with advected values
            Array.Copy(this.newU, this.u, this.u.Length);
            Array.Copy(this.newV, this.v, this.v.Length);
        }

        public void AdvectSmoke(float dt)
        {
            // Copy current smoke density to a temporary array
            this.newM = (double[])this.m.Clone();

            int n = this.numY;
            float h = this.h;
            float h2 = 0.5f * h;

            for (int i = 1; i < this.numX - 1; i++)
            {
                for (int j = 1; j < this.numY - 1; j++)
                {
                    if (this.s[i * n + j] != 0.0f)
                    {
                        double u = (this.u[i * n + j] + this.u[(i + 1) * n + j]) * 0.5f;
                        double v = (this.v[i * n + j] + this.v[i * n + j + 1]) * 0.5f;
                        double x = i * h + h2 - dt * u;
                        double y = j * h + h2 - dt * v;

                        this.newM[i * n + j] = this.SampleField(x, y, SFIELD);

                    }
                }
            }

            // Update smoke density with advected values
            Array.Copy(this.newM, this.m, this.m.Length);
        }

        public void addHorizontalForce()
        {
            for (int i = 0; i < numY; i++)
            {
                if (i < numY * 2 / 5 || i > numY * 3 / 5)
                    continue;
                m[2 * numY + i] = 1;
                u[2 * numY + i] = fluidShower.horizontalForce;
            }
        }

        public void addCircleFluidToMiddle()
        {
            Vector2 Center = new Vector2((int)(numX / 2), (int)(numY / 2));
            float radius = 10f;

            for (int i = 0; i < numX; i++)
            {
                for (int j = 0; j < numY; j++)
                {
                    if ((math.pow(i - Center.x, 2) + math.pow(j - Center.y, 2)) < radius)
                    {
                        m[i * numY + j] = 1;
                    }
                }
            }
        }

        public void setUp()
        {
            for (int i = 0; i < numY; i++)
            {
                if (i < numY * 1 / 3 || i > numY * 2 / 3)
                    continue;
                m[2 * numY + i] = 1;
                //u[2 * numY + i] = fluidShower.horizontalForce;
            }
        }

        public void Simulate(float dt, float gravity)
        {

            if (fluidShower.activateHorizontalFore)
            {
                addHorizontalForce();
            }

            //Array.Clear(p, 0, p.Length);

            AdvectVel(dt);
            AdvectSmoke(dt);

            Integrate(dt, gravity);

            SolveIncompressibility(dt);
            Extrapolate();


        }

        public void setBlocksAsObsticles(HashSet<GameObject> blocks, Vector2 obstacleVelocity)
        {
            /*Given a list of coordinates, puts the s values of the coordinates to 0*/
            int2 coordinate;
            foreach (GameObject block in blocks)
            {
                if (!fluidShower.blockInfos.ContainsKey(block))
                    continue;

                coordinate = fluidShower.blockInfos[block];
                s[coordinate.x * numY + coordinate.y] = 0;
                m[coordinate.x * numY + coordinate.y] = 0;
                u[coordinate.x * numY + coordinate.y] = obstacleVelocity.x;
                v[coordinate.x * numY + coordinate.y] = obstacleVelocity.y;
            }
        }

        public void setBlocksAsNonObsticles(HashSet<GameObject> blocks)
        {
            /*Given a list of coordinates, puts the s values of the coordinates to 0*/
            int2 coordinate;
            foreach (GameObject block in blocks)
            {
                if (!fluidShower.blockInfos.ContainsKey(block))
                    continue;
                Debug.Log("Debolcked object");
                coordinate = fluidShower.blockInfos[block];
                s[coordinate.x * numY + coordinate.y] = 1;
            }
        }

        public void printVolume()
        {
            double totalVolume = 0;

            for(int i = 0; i < numCells; i++)
            {
                totalVolume += m[i];
            }

            Debug.Log(totalVolume);            
        }

        public void printGrid()
        {
            switch (fluidShower.currentField)
            {
                case FieldType.MFIELD:
                    Debug.Log("this.u array: " + string.Join(", ", m));
                    break;
                case FieldType.PFIELD:
                    Debug.Log("this.p array: " + string.Join(", ", p));
                    break;
                case FieldType.UFIELD:
                    Debug.Log("this.u array: " + string.Join(", ", u));
                    break;
                case FieldType.VFIELD:
                    Debug.Log("this.u array: " + string.Join(", ", v));
                    break;
            }
        }
    }


}
