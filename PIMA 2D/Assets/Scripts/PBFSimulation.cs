using System.Collections.Generic;
using UnityEngine;
using Unity.Mathematics;
using UnityEngine.UI;
using JetBrains.Annotations;
using Unity.VisualScripting;
using UnityEditor.UIElements;
using System;

public class PBFSimulation : MonoBehaviour
{
    public GameObject particleObject;

    public int particleNumber;
    public Vector2 particleDistances = new Vector2(1,1);
    public float particleRadius = 1f;
    public float collisionDamping = 0.9f;
    public float gravity = -9.8f;
    public float particleMass = 1;
    public float density = 1;
    public float smoothingRadius = 5f;
    public float constraintForce = 0.01f; //used for stability when we calculate lamda
    public int solverIteration = 1; //the number of times we solve the constraint
    [SerializeField] private SpawnType spawnType;

    private List<GameObject> particlesObjects = new List<GameObject>();

    private Vector2 horizontalBoundaries = new Vector2(0, 10);
    private Vector2 verticalBoundaries = new Vector2(0, 10);

    private Vector3[] positions;
    private Vector3[] predictedPositions;
    private Vector3[] velocities;
    private float[] densities;
    private float[] lamdas;
    private Vector3[] deltaPs;

    private Vector3 gravityVector;

    private Dictionary<int, List<int>> particleNeighbors = new Dictionary<int, List<int>>();

    enum SpawnType { UNIFORM, RANDOM};

    public void Start()
    {
        //particle Data initialization
        positions = new Vector3[particleNumber];
        predictedPositions = new Vector3[particleNumber];
        velocities = new Vector3[particleNumber];
        densities = new float[particleNumber];
        lamdas = new float[particleNumber];
        deltaPs = new Vector3[particleNumber];

        gravityVector = new Vector3(0, gravity, 0);

        //set the boundaries as camera boundaries
        Vector3 cameraCenter = Camera.main.transform.position;
        float halfHeight = Camera.main.orthographicSize;
        float halfWidth = Camera.main.aspect * halfHeight;
        verticalBoundaries = new Vector2(cameraCenter.y - halfHeight, cameraCenter.y + halfHeight);
        horizontalBoundaries = new Vector2(cameraCenter.x - halfWidth, cameraCenter.x + halfWidth);

        if (spawnType == SpawnType.UNIFORM)
            spawnParticles();
        else if (spawnType == SpawnType.RANDOM)
            spawnRandomParticles();

        //initiallize neighbors
        for (int i = 0; i < particleNumber; i++)
        {
            List<int> neighbors = new List<int>();
            particleNeighbors.Add(i, neighbors);
        }
    }

    public void Update()
    {
        //add external forces
        for(int i = 0; i < particleNumber; i++)
        {
            velocities[i] += getExternalForce() * Time.deltaTime;
            predictedPositions[i] = positions[i] + velocities[i] * Time.deltaTime;
        }

        //Find neighbors
        findNeighboors();

        //calculate densities
        for (int i = 0; i < particleNumber; i++)
        {
            densities[i] = calculateDensity(i);
        }

        for(int solveIteration = 0; solveIteration < solverIteration; solveIteration ++ )
        {
            //calculate lamda values
            calculateLamdas();

            //calculate deltaP
            calculateDeltaPs();

            //update predicted pos by respecting the constrain
            for (int i = 0; i < particleNumber; i++)
            {
                predictedPositions[i] += deltaPs[i]*Time.deltaTime;
            }
        }

        //update position and velocity
        for (int i = 0; i < particleNumber; i++)
        {
            velocities[i] = (predictedPositions[i] - positions[i]) / Time.deltaTime;
            positions[i] = predictedPositions[i];
            checkBoundaryConditions(i);
            particlesObjects[i].transform.position = positions[i];
        }

        clearNeighbors();
    }

    private void spawnParticles()
    {

        Vector3 cameraCenter = Camera.main.transform.position;
        Vector3 spawnPos;

        int numColumn = (int)Mathf.Sqrt(particleNumber);
        int lastLineParticles = particleNumber % numColumn;

        Vector3 offset = new Vector3(numColumn / 2 * particleDistances.x, numColumn / 2 * particleDistances.y);

        for (int i = 0; i < particleNumber; i++)
        {
            int row = i / numColumn;
            int column = i % numColumn;

            spawnPos = new Vector3(cameraCenter.x + column * particleDistances.x, cameraCenter.y + row * particleDistances.y, 0) - offset;
            GameObject p1 = Instantiate(particleObject, spawnPos, Quaternion.identity);
            p1.transform.localScale = new Vector3(particleRadius * 2, particleRadius * 2, particleRadius * 2);
            particlesObjects.Add(p1);
            positions[i] = spawnPos;
            velocities[i] = Vector3.zero;
        }

    }

    public void spawnRandomParticles()
    {
        Vector2 randomSpawnPos;
        for (int i = 0; i < particleNumber; i++)
        {

            randomSpawnPos = new Vector3(UnityEngine.Random.Range(horizontalBoundaries.x, horizontalBoundaries.y), UnityEngine.Random.Range(verticalBoundaries.x, verticalBoundaries.y));
            GameObject p1 = Instantiate(particleObject, randomSpawnPos, Quaternion.identity);
            p1.transform.localScale = new Vector3(particleRadius * 2, particleRadius * 2, particleRadius * 2);
            particlesObjects.Add(p1);

            //intialize particle data
            positions[i] = randomSpawnPos;
            velocities[i] = Vector3.zero;
        }
    }

    public void moveParticle()
    {
        for(int i = 0; i < particleNumber; i++)
        {
            positions[i] += velocities[i] * Time.deltaTime;
            particlesObjects[i].transform.position = positions[i];

            checkBoundaryConditions(i);
        }
    }

    private void calculateDeltaPs()
    {
        for(int i = 0; i < particleNumber; i++)
        {
            deltaPs[i] = calculateDeltaP(i);
        }
    }

    private Vector3 calculateDeltaP(int i)
    {
        Vector3 deltaPi = Vector3.zero;
        float dist;
        Vector3 dir;

        foreach(int j in particleNeighbors[i])
        {
            if (i == j)
                continue;

            dist = (predictedPositions[i] - predictedPositions[j]).magnitude;
            dir = dist == 0 ? getRandomDirection() : (predictedPositions[j] - predictedPositions[i]) / dist;

            deltaPi += (lamdas[i] + lamdas[j]) * gradientSpikyKernel(dist, smoothingRadius) * dir;
        }
        return deltaPi;
    }

    private void calculateLamdas()
    {
        for(int i = 0; i < particleNumber; i++)
        {
            lamdas[i] = calculateLamda(i);
        }
    }

    private float calculateLamda(int i)
    {
        float denominator = 0;

        foreach(int k in particleNeighbors[i])
        {
            denominator += math.pow(math.abs(constraintGradient(i, k)), 2);
        }

        return -constraintFunction(i) / (denominator + constraintForce);
    }

    private float constraintFunction(int i)
    {
        return (densities[i] / density) - 1;
    }

    private float constraintGradient(int i, int k)
    {
        //this function is tocalculate the gradient of the constraint function
        //of the particle i by the neighbor k
        float partialDeriv = 0;
        float dist;

        if(k != i)
        {
            dist = (predictedPositions[i] - predictedPositions[k]).magnitude;
            partialDeriv = -gradientSpikyKernel(dist, smoothingRadius);
        }
        else
        {
            foreach(int j in particleNeighbors[i])
            {
                dist = (predictedPositions[i] - predictedPositions[j]).magnitude;
                partialDeriv += gradientSpikyKernel(dist, smoothingRadius);
            }
        }

        return partialDeriv / density;
    }

    private float calculateDensity(int pi)
    {
        float density = 0;
        float dist;

        foreach(int pj in particleNeighbors[pi])
        {
            dist = (predictedPositions[pi] - predictedPositions[pj]).magnitude;  //TO DO dist is supposed to be squared
            density += particleMass * poly6Kernel(dist, smoothingRadius);
        }

        return density;
    }

    private float poly6Kernel(float r, float h)
    {
        if (r >= 0 && r <= h)
        {
            return (315 / (64 * math.PI * math.pow(h, 9))) * math.pow(h * h - r * r, 3);
        }
        return 0;
    }

    private float spikyKernel(float r, float h)
    {
        if (r >= 0 && r <= h)
        {
            return (15 / (math.PI * math.pow(h, 6))) * math.pow(h - r, 3);
        }
        return 0;
    }

    private float gradientSpikyKernel(float r, float h)
    {
        if (r >= 0 && r <= h)
        {
            return -(45 / (math.PI * math.pow(h, 6))) * math.pow(h - r, 2);
        }
        return 0;
    }

    private void checkBoundaryConditions(int pi)
    {
        //vertical conditions
        if ((predictedPositions[pi].y - particleRadius) < verticalBoundaries.x && velocities[pi].y < 0)
        {
            velocities[pi].y = -velocities[pi].y * collisionDamping;
        }
        if ((predictedPositions[pi].y + particleRadius) > verticalBoundaries.y && velocities[pi].y > 0)
        {
            velocities[pi].y = -velocities[pi].y * collisionDamping;
        }
        //horizontal conditions
        if ((predictedPositions[pi].x - particleRadius) < horizontalBoundaries.x && velocities[pi].x < 0)
        {
            velocities[pi].x = -velocities[pi].x * collisionDamping;
        }
        if ((predictedPositions[pi].x + particleRadius) > horizontalBoundaries.y && velocities[pi].x > 0)
        {
            velocities[pi].x = -velocities[pi].x * collisionDamping;
        }
    }

    private void findNeighboors()
    {
        for(int i = 0; i < particleNumber; i++)
        {
            for(int j = 0; j < particleNumber; j++)
            {
                if ((predictedPositions[i] - predictedPositions[j]).magnitude < smoothingRadius)
                {
                    particleNeighbors[i].Add(j);
                }
            }
        }
    }

    private Vector3 getExternalForce()
    {
        return gravityVector;
    }

    private void clearNeighbors()
    {
        for (int i = 0; i < particleNumber; i++)
        {
            List<int> neighbors = new List<int>();
            particleNeighbors[i] = neighbors;
        }
    }

    private Vector3 getRandomDirection()
    {
        return new Vector3(UnityEngine.Random.Range(-1.0f, 1.0f), UnityEngine.Random.Range(-1.0f, 1.0f), 0).normalized;
    }

    private void debugNeighboors(int pi)
    {
        foreach(int i in particleNeighbors[pi])
        {
            particlesObjects[i].GetComponent<SpriteRenderer>().color = Color.red;
        }
    }
}
