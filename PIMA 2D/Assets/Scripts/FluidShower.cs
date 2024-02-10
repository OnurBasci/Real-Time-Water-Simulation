using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using UnityEngine;
using static UnityEngine.Rendering.DebugUI;

public class FluidShower : MonoBehaviour
{
    public GameObject squarePrefab; // Prefab for the grid square
    public float squareSize = 1f; // Size of each square

    public float gravity = -9.81f;
    public float overRelaxation = 1.9f;
    public float density = 1000f;
    public int numX = 100;
    public int numY = 100;
    public float h = 1;
    public float dt = 1.0f/120.0f;
    public int numIter = 100;
    public FieldType currentField;
    public enum FieldType
    {
        MFIELD,
        PFIELD,
        UFIELD,
        VFIELD
    }

    public const int UFIELD = 0;
    public const int VFIELD = 1;
    public const int SFIELD = 2; //smoke field

    private Fluid fluid;

    public void Start()
    {
        fluid = new Fluid(density, numX, numY, h, this);
        GenerateGrid(fluid);
    }

    public void Update()
    {
        fluid.Simulate(dt, gravity, numIter);
        GenerateGrid(fluid);
        fluid.printGrid();
    }


    void GenerateGrid(Fluid fluid)
    {
        float[] gridArray = fluid.m;
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
        float max_pressure = gridArray.Max();
        //if (max_pressure <= 0) max_pressure = gridArray.Min();
        for (int i = 0; i < numX; i++)
        {
            for (int j = 0; j < numY; j++)
            {
                Vector3 spawnPosition = new Vector3(j * squareSize, i * squareSize, 0);
                GameObject square = Instantiate(squarePrefab, spawnPosition, Quaternion.identity);
                square.transform.localScale = new Vector3(squareSize, squareSize, 1f);

                // Set color based on gridArray value
                float c = gridArray[i* numY + j] / max_pressure;

                Color color = new Color(c, c, c);

                square.GetComponent<SpriteRenderer>().color = color;
            }
        }
    }


    public class Fluid
    {
        public FluidShower fluidShower;
        public float density;
        public int numX;
        public int numY;
        public int numCells;
        public float h;
        public float[] u;
        public float[] v;
        public float[] newU;
        public float[] newV;
        public float[] p;
        public float[] s;
        public float[] m;  // m for density (mase)
        public float[] newM;

        public Fluid(float density, int numX, int numY, float h, FluidShower fluidShower)
        {
            this.density = density;
            this.numX = numX;
            this.numY = numY;
            this.numCells = this.numX * this.numY;
            this.h = h;
            this.u = new float[this.numCells];
            this.v = new float[this.numCells];
            this.newU = new float[this.numCells];
            this.newV = new float[this.numCells];
            this.p = new float[this.numCells];
            this.s = new float[this.numCells];
            this.m = new float[this.numCells];
            this.newM = new float[this.numCells];
            this.fluidShower = fluidShower;

            // Initialize m and s array with 1.0
            for (int i = 0; i < numX; i++)
            {
                for(int j = 0; j < numY; j++)
                {
                    //this.m[i] = 1.0f;
                    this.s[i * numY + j] = 1f;
                    if (i == 0 || i == numX - 1 || j == 0)
                        s[i*numY + j] = 0.0f;	// make the border solid
                }
            }

            m[5 * numY + 10] = 1;

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

        private void SolveIncompressibility(int numIters, float dt)
        {
            int n = this.numY;
            float cp = this.density * this.h / dt;

            for (int iter = 0; iter < numIters; iter++)
            {
                for (int i = 1; i < this.numX - 1; i++)
                {
                    for (int j = 1; j < this.numY - 1; j++)
                    {
                        if (this.s[i * n + j] == 0.0f)
                            continue;

                        float s = this.s[i * n + j];
                        float sx0 = this.s[(i - 1) * n + j];
                        float sx1 = this.s[(i + 1) * n + j];
                        float sy0 = this.s[i * n + j - 1];
                        float sy1 = this.s[i * n + j + 1];
                        float sum = sx0 + sx1 + sy0 + sy1;

                        if (sum == 0.0f)
                            continue;

                        float div = this.u[(i + 1) * n + j] - this.u[i * n + j] +
                                    this.v[i * n + j + 1] - this.v[i * n + j];

                        float p = -div / sum;
                        p *= fluidShower.overRelaxation;
                        this.p[i * n + j] += cp * p;

                        this.u[i * n + j] -= sx0 * p;
                        this.u[(i + 1) * n + j] += sx1 * p;
                        this.v[i * n + j] -= sy0 * p;
                        this.v[i * n + j + 1] += sy1 * p;
                    }
                }
            }
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

        public float SampleField(float x, float y, int field)
        {
            //this function applies a biliniar interpolation to calculate the value of a point in a specified field
            int n = this.numY;
            float h = this.h;
            float h1 = 1.0f / h;
            float h2 = 0.5f * h;

            x = Mathf.Max(Mathf.Min(x, this.numX * h), h);
            y = Mathf.Max(Mathf.Min(y, this.numY * h), h);

            float dx = 0.0f;
            float dy = 0.0f;

            float[] f = null;

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

            int x0 = Mathf.Min(Mathf.FloorToInt((x - dx) * h1), this.numX - 1);
            float tx = ((x - dx) - x0 * h) * h1;
            int x1 = Mathf.Min(x0 + 1, this.numX - 1);

            int y0 = Mathf.Min(Mathf.FloorToInt((y - dy) * h1), this.numY - 1);
            float ty = ((y - dy) - y0 * h) * h1;
            int y1 = Mathf.Min(y0 + 1, this.numY - 1);

            float sx = 1.0f - tx;
            float sy = 1.0f - ty;

            float val = sx * sy * f[x0 * n + y0] +
                        tx * sy * f[x1 * n + y0] +
                        tx * ty * f[x1 * n + y1] +
                        sx * ty * f[x0 * n + y1];

            return val;
        }

        public float AvgU(int i, int j)
        {
            int n = this.numY;
            float u = (this.u[i * n + j - 1] + this.u[i * n + j] +
                       this.u[(i + 1) * n + j - 1] + this.u[(i + 1) * n + j]) * 0.25f;
            return u;
        }

        public float AvgV(int i, int j)
        {
            int n = this.numY;
            float v = (this.v[(i - 1) * n + j] + this.v[i * n + j] +
                       this.v[(i - 1) * n + j + 1] + this.v[i * n + j + 1]) * 0.25f;
            return v;
        }

        public void AdvectVel(float dt)
        {
            // Copy current velocities to temporary arrays
            this.newU = (float[])this.u.Clone();
            this.newV = (float[])this.v.Clone();

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
                        float x = i * h;
                        float y = j * h + h2;
                        float u = this.u[i * n + j];
                        float v = this.AvgV(i, j);
                        x = x - dt * u;
                        y = y - dt * v;
                        u = this.SampleField(x, y, UFIELD);
                        this.newU[i * n + j] = u;
                    }
                    // v component
                    if (this.s[i * n + j] != 0.0f && this.s[i * n + j - 1] != 0.0f && i < this.numX - 1)
                    {
                        float x = i * h + h2;
                        float y = j * h;
                        float u = this.AvgU(i, j);
                        float v = this.v[i * n + j];
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
            this.newM = (float[])this.m.Clone();

            int n = this.numY;
            float h = this.h;
            float h2 = 0.5f * h;

            for (int i = 1; i < this.numX - 1; i++)
            {
                for (int j = 1; j < this.numY - 1; j++)
                {
                    if (this.s[i * n + j] != 0.0f)
                    {
                        float u = (this.u[i * n + j] + this.u[(i + 1) * n + j]) * 0.5f;
                        float v = (this.v[i * n + j] + this.v[i * n + j + 1]) * 0.5f;
                        float x = i * h + h2 - dt * u;
                        float y = j * h + h2 - dt * v;

                        this.newM[i * n + j] = this.SampleField(x, y, SFIELD);

                    }
                }
            }

            // Update smoke density with advected values
            Array.Copy(this.newM, this.m, this.m.Length);
        }

        public void addHorizontalForce()
        {
            for(int i = 0; i < numX; i++)
            {
                u[i * numY] = +10;
                //u[i * numY+1] = 10;
                m[i * numY] = 1;
            }
        }

        public void Simulate(float dt, float gravity, int numIters)
        {
            Integrate(dt, gravity);
            //addHorizontalForce();

            Array.Clear(p, 0, p.Length);
            SolveIncompressibility(numIters, dt);

            Extrapolate();
            AdvectVel(dt);
            AdvectSmoke(dt);
            
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
