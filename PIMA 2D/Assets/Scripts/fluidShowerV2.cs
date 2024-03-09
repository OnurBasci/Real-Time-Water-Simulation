using System;
using System.Collections.Generic;
using System.Linq;
using Unity.Mathematics;
using UnityEngine;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;
using static UnityEngine.ParticleSystem;

public class fluidShowerV2 : MonoBehaviour
{
    public GameObject squarePrefab; // Prefab for the grid square
    public GameObject particleObject;
    public GameObject redArrow;
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
    public bool simulateWater = true;

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
    private List<GameObject> particleObjects = new List<GameObject>();
    private List<GameObject> horizontalRedArrows = new List<GameObject>();
    private List<GameObject> verticalRedArrows = new List<GameObject>();

    [NonSerialized] public Fluid fluid;

    const int SOLID = 0, FLUID = 1, EMPTY = 2;

    public void Start()
    {
        blockObjectsArray = new GameObject[numX, numY];

        fluid = new Fluid(density, numX, numY, h, 0.01f, this);
        GenerateGrid(fluid);
        putRedArrowsOnGrid(fluid);
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
        //deltaTime += (Time.deltaTime - deltaTime) * 0.1f;
        fluid.dt = (Time.deltaTime - deltaTime) * 0.1f;

        fluid.Simulate(gravity);
        
        if(simulateWater)
        {
            showWater(fluid);
            updateArrows(fluid);
        }
        else
        {
            simulateGrid(fluid);
        }

        if(printVolume)
        {
            fluid.printVolume();
        }

        if(printGrid)
        {
            fluid.printGrid();
        }
    }

    #region FUNCTION TO PLOT THE GRID

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

    public void showWater(Fluid fluid)
    {
        //this method changes the color of a block by it's type (black for solid, white for air and blue for water)

        for(int i = 0; i < numY; i ++)
        {
            for(int j = 0; j < numX; j++)
            {
                if (fluid.types[j*numY + i] == SOLID)
                {
                    blockObjectsArray[j, i].GetComponent<SpriteRenderer>().color = Color.black;
                }
                else if(fluid.types[j * numY + i] == FLUID)
                {
                    blockObjectsArray[j, i].GetComponent<SpriteRenderer>().color = Color.blue;
                }
                else if(fluid.types[j * numY + i] == EMPTY)
                {
                    blockObjectsArray[j, i].GetComponent<SpriteRenderer>().color = Color.white;
                }
            }
        }
    }

    void putRedArrowsOnGrid(Fluid fluid)
    {
        float h2 = (float)(blockObjectsArray[1, 0].transform.position.x - blockObjectsArray[0, 0].transform.position.x) / (float)2;

        //add horizontal arrows representing horizontal vector field
        for (int i = 0; i < numX + 1; i++)
        {
            for (int j = 0; j < numY; j++)
            {
                if (i == 0)
                {
                    horizontalRedArrows.Add(Instantiate(redArrow, new Vector2(blockObjectsArray[0, j].transform.position.x - h2,
                    blockObjectsArray[0, j].transform.position.y), Quaternion.identity));
                    continue;
                }
                horizontalRedArrows.Add(Instantiate(redArrow, new Vector2(blockObjectsArray[i - 1, j].transform.position.x + h2,
                    blockObjectsArray[i - 1, j].transform.position.y), Quaternion.identity));
            }
        }
        //add vertical arrows representing horizontal vector field
        for (int i = 0; i < numX; i++)
        {
            for (int j = 0; j < numY + 1; j++)
            {
                if (j == 0)
                {
                    verticalRedArrows.Add(Instantiate(redArrow, new Vector2(blockObjectsArray[i, 0].transform.position.x,
                    blockObjectsArray[j, 0].transform.position.y - h2), Quaternion.Euler(0, 0, 90)));
                    continue;
                }
                verticalRedArrows.Add(Instantiate(redArrow, new Vector2(blockObjectsArray[i, j - 1].transform.position.x,
                    blockObjectsArray[i, j - 1].transform.position.y + h2), Quaternion.Euler(0, 0, 90)));
            }
        }
    }

    public void updateArrows(Fluid fluid)
    {
        //update the horizontal velocity field
        float arrowLength = 0;
        float maxVel = FindMax(fluid.u);
        float minVel = FindMin(fluid.u);

        
        for (int i = 0; i < numX; i++)
        {
            for (int j = 0; j < numY; j++)
            {
                if (fluid.u[i * numY + j] > 0)
                {
                    if (maxVel == 0) continue;
                    arrowLength = (float)fluid.u[i * numY + j] / maxVel;
                    horizontalRedArrows[i * numY + j].transform.localScale = new Vector3(arrowLength, 1, 1);
                }
                else if (fluid.u[i * numY + j] < 0)
                {
                    if (minVel == 0) continue;
                    arrowLength = (float)fluid.u[i * numY + j] / minVel;
                    horizontalRedArrows[i * numY + j].transform.localScale = new Vector3(-arrowLength, 1, 1);
                }
            }
        }

        //update the vertical field
        maxVel = FindMax(fluid.v);
        minVel = FindMin(fluid.v);
        for (int i = 0; i < numX; i++)
        {
            for (int j = 0; j < numY; j++)
            {
                if (fluid.v[i * numY + j] > 0)
                {
                    if (maxVel == 0) continue;
                    arrowLength = (float)fluid.v[i*numY + j] / maxVel;
                    verticalRedArrows[i *numY + j].transform.localScale = new Vector3(arrowLength, 1, 1);
                }
                else if (fluid.v[i * numY + j] < 0)
                {
                    if (minVel == 0) continue;
                    arrowLength = (float)fluid.v[i * numY + j] / minVel;
                    verticalRedArrows[i * numY + j].transform.localScale = new Vector3(-arrowLength, 1, 1);
                }
            }
        }
    }

    #endregion

    #region USEFULL FUNCTION

    public static int FindMax(double[] array)
    {
        int max = int.MinValue;
        foreach (int element in array)
        {
            if (element > max)
                max = element;
        }
        return max;
    }

    public static int FindMin(double[] array)
    {
        int min = int.MaxValue;
        foreach (int element in array)
        {
            if (element < min)
                min = element;
        }
        return min;
    }

    #endregion

    public class Fluid
    {
        public fluidShowerV2 fluidShower;
        public float density;
        public float dt;
        public int numX;
        public int numY;
        public int numCells;
        public float h;
        public double[] u;
        public double[] v;
        public double[] newU;
        public double[] newV;
        public double[] p;
        public int[] types; //takes 0 for solid 1 for fluid and 2 for air
        public double[] m;  // m for density (mase)
        public double[] newM;

        List<Vector2> particles;

        public Fluid(float density, int numX, int numY, float h, float dt, fluidShowerV2 fluidShower)
        {
            this.density = density;
            this.dt = dt;
            this.numX = numX;
            this.numY = numY;
            this.numCells = this.numX * this.numY;
            this.h = h;
            this.u = new double[this.numCells];
            this.v = new double[this.numCells];
            this.newU = new double[this.numCells];
            this.newV = new double[this.numCells];
            this.p = new double[this.numCells];
            this.types = new int[this.numCells];
            this.m = new double[this.numCells];
            this.newM = new double[this.numCells];
            this.fluidShower = fluidShower;

            // Initialize types
            for (int i = 0; i < numX; i++)
            {
                for (int j = 0; j < numY; j++)
                {
                    types[i * numY + j] = EMPTY;
                }
            }

            for (int i = 0; i < numX; i++)
            {
                types[0 * numY + i] = SOLID;
                types[i * numY + 0] = SOLID;
                types[(numX - 1) *numY + i] = SOLID;
                types[i * numY + (numY - 1)] = SOLID;
            }

            for (int j = 1; j < numY - 1; j++)
            {
                for (int i = 1; i < numX / 2; i++)
                {
                    types[i* numY + j] = FLUID;
                }
            }

            //initialize particles
            this.particles = new List<Vector2>();
            for (int i = 0; i < numX; i++)
            {
                for (int j = 0; j < numY; j++)
                {
                    if (types[i * numY + j] == FLUID)
                    {
                        Vector2 randomOffset1 = new Vector2(UnityEngine.Random.Range(-0.25f, 0.25f), UnityEngine.Random.Range(-0.25f, 0.25f));
                        Vector2 randomOffset2 = new Vector2(UnityEngine.Random.Range(-0.25f, 0.25f), UnityEngine.Random.Range(-0.25f, 0.25f));
                        Vector2 randomOffset3 = new Vector2(UnityEngine.Random.Range(-0.25f, 0.25f), UnityEngine.Random.Range(-0.25f, 0.25f));
                        Vector2 randomOffset4 = new Vector2(UnityEngine.Random.Range(-0.25f, 0.25f), UnityEngine.Random.Range(-0.25f, 0.25f));
                        particles.Add(new Vector2(i - 0.25f, j - 0.25f) + randomOffset1);
                        particles.Add(new Vector2(i - 0.25f, j - 0.75f) + randomOffset2);
                        particles.Add(new Vector2(i + 0.25f, j + 0.25f) + randomOffset3);
                        particles.Add(new Vector2(i + 0.25f, j + 0.25f) + randomOffset4);
                    }
                }
            }

            foreach (var particle in particles)
            {
                var pObject = Instantiate(fluidShower.particleObject, particle, Quaternion.identity);
                fluidShower.particleObjects.Add(pObject);
            }

        }

        private void Integrate(float gravity)
        {
            int n = this.numY;
            for (int i = 1; i < this.numX; i++)
            {
                for (int j = 1; j < this.numY - 1; j++)
                {
                    if (this.types[i * n + j] != 0.0f && this.types[i * n + j - 1] != 0.0f)
                    {
                        this.v[i * n + j] += gravity * dt;
                    }
                }
            }
        }

        private void SolveIncompressibility()
        {
            DenseMatrix u = DenseMatrix.Create(numX, numY, 0);
            DenseMatrix v = DenseMatrix.Create(numX, numY, 0);
            DenseMatrix p = DenseMatrix.Create(numX, numY, 0);
            DenseMatrix t = DenseMatrix.Create(numX, numY, 0);

            for(int i = 0; i < numX; i ++)
            {
                for(int j = 0; j < numY; j++)
                {
                    u[i,j] = this.u[i*numY + j];
                    v[i,j] = this.v[i*numY + j];
                    p[i,j] = this.p[i*numY + j];
                    t[i,j] = this.types[i*numY + j];
                }
            }

            Projector projector = new Projector(u, v, p, t, numX, dt, density);
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


        void Project()
        {
            List<Vector2Int> fluidCells = new List<Vector2Int>();
            Dictionary<Vector2Int, int> fluidDict = new Dictionary<Vector2Int, int>();
            for (int i = 0; i < numX; i++)
            {
                for (int j = 0; j < numY; j++)
                {
                    if (types[i * numY + j] == FLUID)
                    {
                        fluidCells.Add(new Vector2Int(i, j));
                        fluidDict[new Vector2Int(i, j)] = fluidCells.Count - 1;
                    }
                }
            }
            int nFluid = fluidCells.Count;
            //float[,] A = new float[nFluid, nFluid];
            //float[] b = new float[nFluid];

            var A = Matrix<float>.Build.Dense(nFluid, nFluid);
            var b = Vector<float>.Build.Dense(nFluid);

            for (int idx = 0; idx < nFluid; idx++)
            {
                int i = fluidCells[idx].x;
                int j = fluidCells[idx].y;

                if (types[(i - 1) * numY + j] == FLUID)
                    A[idx, fluidDict[new Vector2Int(i - 1, j)]] = -1;
                if (types[i * numY + j - 1] == FLUID)
                    A[idx, fluidDict[new Vector2Int(i, j - 1)]] = -1;
                if (types[(i + 1) * numY + j] == FLUID)
                    A[idx, fluidDict[new Vector2Int(i + 1, j)]] = -1;
                if (types[i * numY + j + 1] == FLUID)
                    A[idx, fluidDict[new Vector2Int(i, j + 1)]] = -1;

                int nonSolid = 0;
                double D = 0;

                if (types[(i - 1) * numY + j] != SOLID)
                {
                    D -= u[i * numY + j];
                    nonSolid++;
                }
                if (types[i * numY + j - 1] != SOLID)
                {
                    D -= v[i * numY + j];
                    nonSolid++;
                }
                if (types[(i + 1) * numY + j] != SOLID)
                {
                    D += u[(i + 1) * numY + j];
                    nonSolid++;
                }
                if (types[i * numY + j + 1] != SOLID)
                {
                    D += v[i * numY + j + 1];
                    nonSolid++;
                }

                A[idx, idx] = nonSolid;
                b[idx] = -(density / dt) * (float)D;
            }

            //float[] x = new float[nFluid];
            // Use a solver to solve Ax = b, you can use a solver library or implement your own solver
            var x = A.Solve(b);

            // Update u2
            double[] u2 = new double[u.Length];
            Array.Copy(u, u2, u.Length);

            for (int i = 0; i < numX - 1; i++)
            {
                for (int j = 0; j < numY; j++)
                {
                    if (types[i * numY + j] != FLUID && types[(i + 1) * numY + j] != FLUID)
                        continue;

                    float rightP = 0, leftP = 0;

                    if (types[(i + 1) * numY + j] == FLUID)
                        rightP = x[fluidDict[new Vector2Int(i + 1, j)]];
                    if (types[(i + 1) * numY + j] == EMPTY)
                        rightP = 0;
                    if (types[(i + 1) * numY + j] == SOLID)
                    {
                        u2[(i + 1) * numY + j] = 0;
                        continue;
                    }

                    if (types[i * numY + j] == FLUID)
                        leftP = x[fluidDict[new Vector2Int(i, j)]];
                    if (types[i * numY + j] == EMPTY)
                        leftP = 0;
                    if (types[i * numY + j] == SOLID)
                    {
                        u2[(i + 1) * numY + j] = 0;
                        continue;
                    }

                    u2[(i + 1) * numY + j] -= (dt / density * (rightP - leftP));
                }
            }

            // Update v2
            double[] v2 = new double[v.Length];
            Array.Copy(v, v2, v.Length);

            for (int i = 0; i < numX; i++)
            {
                for (int j = 0; j < numX - 1; j++)
                {
                    if (types[i * numY + j] != FLUID && types[i * numY + j + 1] != FLUID)
                        continue;

                    float topP = 0, bottomP = 0;

                    if (types[i * numY + j + 1] == FLUID)
                        topP = x[fluidDict[new Vector2Int(i, j + 1)]];
                    if (types[i * numY + j + 1] == EMPTY)
                        topP = 0;
                    if (types[i * numY + j + 1] == SOLID)
                    {
                        v2[i * numY + j + 1] = 0;
                        continue;
                    }

                    if (types[i * numY + j] == FLUID)
                        bottomP = x[fluidDict[new Vector2Int(i, j)]];
                    if (types[i * numY + j] == EMPTY)
                        bottomP = 0;
                    if (types[i * numY + j] == SOLID)
                    {
                        v2[i * numY + j + 1] = 0;
                        continue;
                    }

                    v2[i * numY + j + 1] -= (dt / density * (topP - bottomP));
                }
            }
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

        public void AdvectVel()
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
                    if (this.types[i * n + j] != 0.0f && this.types[(i - 1) * n + j] != 0.0f && j < this.numY - 1)
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
                    if (this.types[i * n + j] != 0.0f && this.types[i * n + j - 1] != 0.0f && i < this.numX - 1)
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

        public void AdvectSmoke()
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
                    if (this.types[i * n + j] != 0.0f)
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

        void MoveParticles()
        {
            // Clear the types array where it's not SOLID
            for (int i = 0; i < numX; i++)
            {
                for (int j = 0; j < numY; j++)
                {
                    if (types[i * numY + j] != SOLID)
                    {
                        types[i * numY + j] = EMPTY;
                    }
                }
            }

            // Move particles using RK2 method
            for (int idx = 0; idx < particles.Count; idx++)
            {
                Vector2 newParticlePosition = RK2(particles[idx]);
                float newX = newParticlePosition.x;
                float newY = newParticlePosition.y;

                // Check if the new position is not SOLID, then update particles and types array
                if (types[(int)newX * numY + (int)newY] != SOLID)
                {
                    //Debug.Log("before at " + (newX, newY) + " " + types[newX, newY]);
                    particles[idx] = newParticlePosition;
                    types[(int)newX * numY + (int)newY] = FLUID;

                    //Debug.Log("after at " + (newX, newY) + " " + types[newX, newY]);
                }
            }

            for (int idx = 0; idx < particles.Count; idx++)
            {
                fluidShower.particleObjects[idx].transform.position = particles[idx];
            }
        }

        float BilinearInterpU(Vector2 p)
        {
            float x = Mathf.Clamp(p.x, 0, numX - 0.0001f);
            float y = Mathf.Clamp(p.y, 0.5f, numY - 0.5001f);
            int x1 = (int)x;
            int x2 = x1 + 1;
            int y1 = (int)(y - 0.5f) + 1;
            int y2 = y1 + 1;
            int i = (int)x;
            int j = (int)(y - 0.5f);

            var b = Matrix<float>.Build.Dense(2, 2);

            b[0, 0] = (float)u[i * numY + j];
            b[0, 1] = (float)u[i * numY + j + 1];
            b[1, 0] = (float)u[(i + 1) * numY + j];
            b[1, 1] = (float)u[(i + 1) * numY + j + 1];

            Vector<float> a = Vector<float>.Build.Dense(new float[] { x2 - x, x - x1 });
            Vector<float> c = Vector<float>.Build.Dense(new float[] { y2 - y, y - y1 });

            return (a.ToRowMatrix() * b * c.ToColumnMatrix())[0, 0] / ((x2 - x1) * (y2 - y1));

            /*
            float Q11 = u[i, j];
            float Q12 = u[i, j + 1];
            float Q21 = u[i + 1, j];
            float Q22 = u[i + 1, j + 1];
            float a0 = x2 - x;
            float a1 = x - x1;
            float b0 = y2 - y;
            float b1 = y - y1;
            return a0 * (b0 * Q11 + b1 * Q12) + a1 * (b0 * Q21 + b1 * Q22);
             */
        }

        float BilinearInterpV(Vector2 p)
        {
            float x = Mathf.Clamp(p.x, 0.5f, numX - 0.5001f);
            float y = Mathf.Clamp(p.y, 0, numY - 0.0001f);
            int x1 = (int)(x - 0.5f) + 1;
            int x2 = x1 + 1;
            int y1 = (int)y;
            int y2 = y1 + 1;
            int i = (int)(x - 0.5f);
            int j = (int)y;


            var b = Matrix<float>.Build.Dense(2, 2);

            b[0, 0] = (float)v[i * numY + j];
            b[0, 1] = (float)v[i * numY + j + 1];
            b[1, 0] = (float)v[(i + 1) * numY + j];
            b[1, 1] = (float)v[(i + 1) * numY + j + 1];

            Vector<float> a = Vector<float>.Build.Dense(new float[] { x2 - x, x - x1 });
            Vector<float> c = Vector<float>.Build.Dense(new float[] { y2 - y, y - y1 });

            return (a.ToRowMatrix() * b * c.ToColumnMatrix())[0, 0] / ((x2 - x1) * (y2 - y1));

            /*float Q11 = v[i, j];
            float Q12 = v[i, j + 1];
            float Q21 = v[i + 1, j];
            float Q22 = v[i + 1, j + 1];
            float a0 = x2 - x;
            float a1 = x - x1;
            float b0 = y2 - y;
            float b1 = y - y1;
            y = a0 * (b0 * Q11 + b1 * Q12) + a1 * (b0 * Q21 + b1 * Q22);
            return y;*/
        }


        Vector2 RK2(Vector2 curr)
        {
            Vector2 k1 = new Vector2(BilinearInterpU(curr), BilinearInterpV(curr));
            Vector2 k2 = new Vector2(BilinearInterpU(curr + 0.5f * dt * k1), BilinearInterpV(curr + 0.5f * dt * k1));
            return curr + dt * k2;
        }

        Vector2 BackwardsRK2(Vector2 curr)
        {
            Vector2 k1 = new Vector2(BilinearInterpU(curr), BilinearInterpV(curr));
            Vector2 k2 = new Vector2(BilinearInterpU(curr - 0.5f * dt * k1), BilinearInterpV(curr - 0.5f * dt * k1));
            return curr - dt * k2;
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



        public void Simulate(float gravity)
        {
            if (fluidShower.activateHorizontalFore)
            {
                addHorizontalForce();
            }

            //Array.Clear(p, 0, p.Length);

            AdvectVel();
            //AdvectSmoke();
            //SolveIncompressibility();
            Project();
            Extrapolate();

            MoveParticles();

            Integrate(gravity);


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
                types[coordinate.x * numY + coordinate.y] = 0;
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
                types[coordinate.x * numY + coordinate.y] = 1;
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
