using MathNet.Numerics.LinearAlgebra;
using System;
using System.Collections;
using System.Collections.Generic;
using System.Drawing;
using System.Linq;
using Unity.Mathematics;
using UnityEngine;
using static UnityEngine.ParticleSystem;

public class FluidSimulationV3 : MonoBehaviour
{
    public GameObject squarePrefab; // Prefab for the grid square
    public GameObject particleObject;
    public GameObject redArrow;
    public float squareSize = 1f; // Size of each square

    public float gravity = -0.4f;
    public float density = 1;
    public int N = 16;
    public float dt = 0.2f;


    float[,] u;
    float[,] v;
    const int SOLID = 0, FLUID = 1, EMPTY = 2;
    int[,] types;

    List<Vector2> particles = new List<Vector2>();

    private float deltaTime;
    private GameObject[,] blockObjectsArray;
    private List<GameObject> particleObjects = new List<GameObject>();
    private List<GameObject> horizontalRedArrows = new List<GameObject>();
    private List<GameObject> verticalRedArrows = new List<GameObject>();

    void Start()
    {

        blockObjectsArray = new GameObject[N, N];

        u = new float[N + 1, N];
        v = new float[N, N + 1];
        types = new int[N, N];

        GenerateGrid();
        putRedArrowsOnGrid();

        InitializeTypes();
        InitializeParticles();

        Project();
    }

    void Update()
    {
        //deltaTime += (Time.deltaTime - deltaTime) * 0.1f;
        //dt = deltaTime;


        Advect();
        Project();
        Extrapolate();
        MoveParticles();

        //LogArray2D(types);

        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N - 1; j++)
            {
                if (types[i, j] == FLUID || types[i, j + 1] == FLUID)
                {
                    v[i, j + 1] += dt * gravity;
                }
            }
        }

        updateGrid();
    }

    void InitializeTypes()
    {
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                types[i, j] = EMPTY;
            }
        }

        for (int i = 0; i < N; i++)
        {
            types[0, i] = SOLID;
            types[i, 0] = SOLID;
            types[N - 1, i] = SOLID;
            types[i, N - 1] = SOLID;
        }

        for (int j = 1; j < N - 1; j++)
        {
            for (int i = 1; i < N / 2; i++)
            {
                types[i, j] = FLUID;
            }
        }
    }

    void InitializeParticles()
    {
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                if (types[i, j] == FLUID)
                {
                    Vector2 randomOffset1 = new Vector2(UnityEngine.Random.Range(-0.25f, 0.25f), UnityEngine.Random.Range(-0.25f, 0.25f));
                    Vector2 randomOffset2 = new Vector2(UnityEngine.Random.Range(-0.25f, 0.25f), UnityEngine.Random.Range(-0.25f, 0.25f));
                    Vector2 randomOffset3 = new Vector2(UnityEngine.Random.Range(-0.25f, 0.25f), UnityEngine.Random.Range(-0.25f, 0.25f));
                    Vector2 randomOffset4 = new Vector2(UnityEngine.Random.Range(-0.25f, 0.25f), UnityEngine.Random.Range(-0.25f, 0.25f));
                    particles.Add(new Vector2(i + 0.25f, j + 0.25f) + randomOffset1);
                    particles.Add(new Vector2(i + 0.25f, j + 0.75f) + randomOffset2);
                    particles.Add(new Vector2(i + 0.75f, j + 0.25f) + randomOffset3);
                    particles.Add(new Vector2(i + 0.75f, j + 0.75f) + randomOffset4);
                }
            }
        }

        foreach (var particle in particles)
        {
            var pObject = Instantiate(particleObject, particle, Quaternion.identity);
            particleObjects.Add(pObject);
        }
    }

    float BilinearInterpU(Vector2 p)
    {
        float x = Mathf.Clamp(p.x, 0, N - 0.0001f);
        float y = Mathf.Clamp(p.y, 0.5f, N - 0.5001f);
        int x1 = (int)x;
        int x2 = x1 + 1;
        int y1 = (int)(y - 0.5f) + 1;
        int y2 = y1 + 1;
        int i = (int)x;
        int j = (int)(y - 0.5f);

        var b = Matrix<float>.Build.Dense(2, 2);

        b[0, 0] = u[i, j];
        b[0, 1] = u[i, j+1];
        b[1, 0] = u[i+1, j];
        b[1, 1] = u[i+1, j+1];

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
        float x = Mathf.Clamp(p.x, 0.5f, N - 0.5001f);
        float y = Mathf.Clamp(p.y, 0, N - 0.0001f);
        int x1 = (int)(x - 0.5f) + 1;
        int x2 = x1 + 1;
        int y1 = (int)y;
        int y2 = y1 + 1;
        int i = (int)(x - 0.5f);
        int j = (int)y;


        var b = Matrix<float>.Build.Dense(2, 2);

        b[0, 0] = v[i, j];
        b[0, 1] = v[i, j + 1];
        b[1, 0] = v[i + 1, j];
        b[1, 1] = v[i + 1, j + 1];

        Vector<float> a = Vector<float>.Build.Dense(new float[] { x2 - x, x - x1 });
        Vector<float> c = Vector<float>.Build.Dense(new float[] { y2 - y, y - y1 });

        return (a.ToRowMatrix() * b * c.ToColumnMatrix())[0, 0]/((x2 - x1)*(y2 - y1));

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

    void Extrapolate()
    {
        // Define maximum integer value for marking unreachable cells
        int maxint = N * N;

        // Array to store distances to nearest FLUID cells
        int[,] d = new int[N + 1, N];
        for(int i = 0; i < N+1; i++)
        {
            for(int j = 0; j < N; j++)
            {
                d[i, j] = maxint;
            }
        }

        // List to store cells to be processed
        List<Tuple<int, int>> W = new List<Tuple<int, int>>();

        // First pass for u
        //put the distance of the inside fluids to 0
        for (int i = 0; i < N - 1; i++)
        {
            for (int j = 0; j < N; j++)
            {
                if (types[i, j] == FLUID || types[i + 1, j] == FLUID)
                {
                    d[i + 1, j] = 0;
                }
            }
        }

        // for each of the blocks that are not fluid if there is a fluid around them put distance to 1
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
                        newU += u[i2, j2];
                        count++;
                    }
                }
            }
            u[i, j] = (float)newU / count;

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
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N + 1; j++)
            {
                d[i, j] = maxint;
            }
        }
        W.Clear();

        // Second pass for v
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N - 1; j++)
            {
                if (types[i, j] == FLUID || types[i, j + 1] == FLUID)
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
                        newV += v[i2, j2];
                        count++;
                    }
                }
            }
            v[i, j] = (float)newV / count;

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

    void Advect()
    {
        float[,] u2 = new float[N + 1, N];
        float[,] v2 = new float[N, N + 1];

        // Advecting horizontal velocity u
        for (int i = 0; i < N + 1; i++)
        {
            for (int j = 0; j < N; j++)
            {
                Vector2 curr = new Vector2(i, j + 0.5f);
                Vector2 prev = BackwardsRK2(curr);
                u2[i, j] = BilinearInterpU(prev);
            }
        }

        // Advecting vertical velocity v
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N + 1; j++)
            {
                Vector2 curr = new Vector2(i + 0.5f, j);
                Vector2 prev = BackwardsRK2(curr);
                v2[i, j] = BilinearInterpV(prev);;
            }
        }

        u = u2;
        v = v2;
    }

    void MoveParticles()
    {
        // Clear the types array where it's not SOLID
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                if (types[i, j] != SOLID)
                {
                    types[i, j] = EMPTY;
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
            if (types[(int)newX, (int)newY] != SOLID)
            {
                //Debug.Log("before at " + (newX, newY) + " " + types[newX, newY]);
                particles[idx] = newParticlePosition;
                types[(int)newX, (int)newY] = FLUID;

                //Debug.Log("after at " + (newX, newY) + " " + types[newX, newY]);
            }
        }

        for (int idx = 0; idx < particles.Count; idx++)
        {
            particleObjects[idx].transform.position = particles[idx];
        }
    }

    void Project()
    {
        List<Vector2Int> fluidCells = new List<Vector2Int>();
        Dictionary<Vector2Int, int> fluidDict = new Dictionary<Vector2Int, int>();
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                if (types[i, j] == FLUID)
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

            if (types[i - 1, j] == FLUID)
                A[idx, fluidDict[new Vector2Int(i - 1, j)]] = -1;
            if (types[i, j - 1] == FLUID)
                A[idx, fluidDict[new Vector2Int(i, j - 1)]] = -1;
            if (types[i + 1, j] == FLUID)
                A[idx, fluidDict[new Vector2Int(i + 1, j)]] = -1;
            if (types[i, j + 1] == FLUID)
                A[idx, fluidDict[new Vector2Int(i, j + 1)]] = -1;

            int nonSolid = 0;
            float D = 0;

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
            b[idx] = -(density / dt) * D;
        }

        //float[] x = new float[nFluid];
        // Use a solver to solve Ax = b, you can use a solver library or implement your own solver
        var x = A.Solve(b);

        // Update u2
        float[,] u2 = new float[N+1, N];
        Array.Copy(u, u2, u.Length);

        for (int i = 0; i < N - 1; i++)
        {
            for (int j = 0; j < N; j++)
            {
                if (types[i, j] != FLUID && types[i + 1, j] != FLUID)
                    continue;

                float rightP = 0, leftP = 0;

                if (types[i + 1, j] == FLUID)
                    rightP = x[fluidDict[new Vector2Int(i + 1, j)]];
                if (types[i + 1, j] == EMPTY)
                    rightP = 0;
                if (types[i + 1, j] == SOLID)
                {
                    u2[i + 1, j] = 0;
                    continue;
                }

                if (types[i, j] == FLUID)
                    leftP = x[fluidDict[new Vector2Int(i, j)]];
                if (types[i, j] == EMPTY)
                    leftP = 0;
                if (types[i, j] == SOLID)
                {
                    u2[i + 1, j] = 0;
                    continue;
                }

                u2[i + 1, j] -= (dt/density * (rightP - leftP));
            }
        }

        // Update v2
        float[,] v2 = new float[N, N+1];
        Array.Copy(v, v2, v.Length);

        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N - 1; j++)
            {
                if (types[i, j] != FLUID && types[i, j + 1] != FLUID)
                    continue;

                float topP = 0, bottomP = 0;

                if (types[i, j + 1] == FLUID)
                    topP = x[fluidDict[new Vector2Int(i, j + 1)]];
                if (types[i, j + 1] == EMPTY)
                    topP = 0;
                if (types[i, j + 1] == SOLID)
                {
                    v2[i, j + 1] = 0;
                    continue;
                }

                if (types[i, j] == FLUID)
                    bottomP = x[fluidDict[new Vector2Int(i, j)]];
                if (types[i, j] == EMPTY)
                    bottomP = 0;
                if (types[i, j] == SOLID)
                {
                    v2[i, j + 1] = 0;
                    continue;
                }

                v2[i, j + 1] -= (dt/density * (topP - bottomP));
            }
        }
        u = u2;
        v = v2;
    }

    void GenerateGrid()
    {

        //if (max_pressure <= 0) max_pressure = gridArray.Min();
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                Vector3 spawnPosition = new Vector3(j * squareSize, i * squareSize, 0);
                GameObject square = Instantiate(squarePrefab, spawnPosition, Quaternion.identity);
                square.transform.localScale = new Vector3(squareSize, squareSize, 1f);

                blockObjectsArray[j, i] = square;
            }
        }
    }

    void putRedArrowsOnGrid()
    {
        float h2 = (float)(blockObjectsArray[1, 0].transform.position.x - blockObjectsArray[0, 0].transform.position.x)/(float)2;

        //add horizontal arrows representing horizontal vector field
        for (int i = 0; i < u.GetLength(0); i++)
        {
            for (int j = 0;  j < u.GetLength(1); j++)
            {
                if (i == 0)
                {
                    horizontalRedArrows.Add(Instantiate(redArrow, new Vector2(blockObjectsArray[0, j].transform.position.x - h2,
                    blockObjectsArray[0, j].transform.position.y), Quaternion.identity));
                    continue;
                }
                horizontalRedArrows.Add(Instantiate(redArrow, new Vector2(blockObjectsArray[i-1, j].transform.position.x + h2, 
                    blockObjectsArray[i-1, j].transform.position.y),Quaternion.identity));
            }
        }
        //add vertical arrows representing horizontal vector field
        for (int i = 0; i < v.GetLength(0); i++)
        {
            for (int j = 0; j < v.GetLength(1); j++)
            {
                if (j == 0)
                {
                    verticalRedArrows.Add(Instantiate(redArrow, new Vector2(blockObjectsArray[i, 0].transform.position.x ,
                    blockObjectsArray[j, 0].transform.position.y - h2), Quaternion.Euler(0,0,90)));
                    continue;
                }
                verticalRedArrows.Add(Instantiate(redArrow, new Vector2(blockObjectsArray[i, j-1].transform.position.x,
                    blockObjectsArray[i, j-1].transform.position.y + h2), Quaternion.Euler(0,0,90)));
            }
        }
    }


    public static void LogArray2D<T>(T[,] array2D)
    {
        int numRows = array2D.GetLength(0);
        int numCols = array2D.GetLength(1);
        string arrayString = "";

        for (int i = 0; i < numRows; i++)
        {
            for (int j = 0; j < numCols; j++)
            {
                arrayString += array2D[i, j];
                if (j < numCols - 1)
                {
                    arrayString += ",";
                }
            }
            arrayString += "\n";
            if (i < numRows - 1)
            {
                arrayString += "; ";
            }
        }

        Debug.Log(arrayString);
    }

    public void updateGrid()
    {
        //update the blocks color
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                if (types[i,j] == 1)
                {
                    blockObjectsArray[i, j].GetComponent<SpriteRenderer>().color = UnityEngine.Color.blue;
                }
                if (types[i, j] == 2)
                {
                    blockObjectsArray[i, j].GetComponent<SpriteRenderer>().color = UnityEngine.Color.white;
                }
                if (types[i, j] == 0)
                {
                    blockObjectsArray[i, j].GetComponent<SpriteRenderer>().color = UnityEngine.Color.black;
                }
            }
        }

        //update the horizontal velocity field
        float arrowLength = 0;
        float maxVel = FindMax(u);
        float minVel = FindMin(u);
        for(int i = 0; i < u.GetLength(0); i++)
        {
            for(int j = 0; j < u.GetLength(1); j++)
            {
                if (u[i,j] > 0)
                {
                    if (maxVel == 0) continue;
                    arrowLength = u[i,j]/maxVel;
                    horizontalRedArrows[i * u.GetLength(1) + j].transform.localScale = new Vector3(arrowLength, 1,1);
                }
                else if(u[i, j] < 0)
                {
                    if (minVel == 0) continue;
                    arrowLength = u[i, j] / minVel;
                    horizontalRedArrows[i * u.GetLength(1) + j].transform.localScale = new Vector3(-arrowLength, 1, 1);
                }
            }
        }

        //update the vertical field
        maxVel = FindMax(u);
        minVel = FindMin(u);
        for (int i = 0; i < v.GetLength(0); i++)
        {
            for (int j = 0; j < v.GetLength(1); j++)
            {
                if (v[i, j] > 0)
                {
                    if (maxVel == 0) continue;
                    arrowLength = v[i, j] / maxVel;
                    verticalRedArrows[i * v.GetLength(1) + j].transform.localScale = new Vector3(arrowLength, 1, 1);
                }
                else if (v[i, j] < 0)
                {
                    if (minVel == 0) continue;
                    arrowLength = v[i, j] / minVel;
                    verticalRedArrows[i * v.GetLength(1) + j].transform.localScale = new Vector3(-arrowLength, 1, 1);
                }
            }
        }

    }

    void printDivergence()
    {
        List<Vector2Int> fluidCells = new List<Vector2Int>();
        Dictionary<Vector2Int, int> fluidDict = new Dictionary<Vector2Int, int>();
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                if (types[i, j] == FLUID)
                {
                    fluidCells.Add(new Vector2Int(i, j));
                    fluidDict[new Vector2Int(i, j)] = fluidCells.Count - 1;
                }
            }
        }
        int nFluid = fluidCells.Count;

        //check if the divergenc is 0 
        List<float> divs = new List<float>();
        for (int idx = 0; idx < nFluid; idx++)
        {
            int i = fluidCells[idx].x;
            int j = fluidCells[idx].y;

            int nonSolid = 0;
            float D = 0;

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
            divs.Add(D);
        }

        float m = float.MinValue;

        foreach (float div in divs)
        {
            if (div > m)
            {
                m = div;
            }
        }
        Debug.Log("max div is " + m);
    }

    public static int FindMax(float[,] array)
    {
        int max = int.MinValue;
        foreach (int element in array)
        {
            if (element > max)
                max = element;
        }
        return max;
    }

    public static int FindMin(float[,] array)
    {
        int min = int.MaxValue;
        foreach (int element in array)
        {
            if (element < min)
                min = element;
        }
        return min;
    }
}
