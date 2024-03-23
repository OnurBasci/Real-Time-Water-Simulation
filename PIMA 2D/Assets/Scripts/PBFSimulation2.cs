using System.Collections.Generic;
using UnityEngine;
using Unity.Mathematics;
using Unity.Jobs;
using Unity.Collections;
using Unity.VisualScripting;
using UnityEngine.Animations;
using MathNet.Numerics.Statistics.Mcmc;

public class PBFSimulation2 : MonoBehaviour
{
    //this class is for position based Simulation. Basicly it is the same as PBFSimulation but it is optimized with Unitys job system

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
    public float deltaTime = 0.016f;
    [SerializeField] private SpawnType spawnType;
    public bool FindNeighborsWithGrid = true;
    //tension stability coefficients
    public float k = 0.1f; //a small positif coefficient
    public float deltaq = 0.1f;
    public bool activateSurfaceTension = false;

    private List<GameObject> particlesObjects = new List<GameObject>();

    private Vector2 horizontalBoundaries = new Vector2(0, 10);
    private Vector2 verticalBoundaries = new Vector2(0, 10);

    private NativeArray<Vector3> positions;
    private NativeArray<Vector3> predictedPositions;
    private NativeArray<Vector3> velocities;
    private NativeArray<float> densities;
    private NativeArray<float> lamdas;
    private NativeArray<Vector3> deltaPs;

    private Vector3 gravityVector;

    private NativeArray<int> particleNeighbors; //this is a particleNumber * particleNumber array to acces the jth neighbor of ith particle particleNeighbor[i*n+j]
    private NativeArray<int> neighborCounter;

    private GridManager gridManager;

    enum SpawnType { UNIFORM, RANDOM};

    public void Start()
    {
        //particle Data initialization
        positions = new NativeArray<Vector3>(particleNumber, Allocator.Persistent);
        predictedPositions = new NativeArray<Vector3>(particleNumber, Allocator.Persistent);
        velocities = new NativeArray<Vector3>(particleNumber, Allocator.Persistent);
        densities = new NativeArray<float>(particleNumber, Allocator.Persistent);
        lamdas = new NativeArray<float>(particleNumber, Allocator.Persistent);
        deltaPs = new NativeArray<Vector3>(particleNumber, Allocator.Persistent);

        gravityVector = new Vector3(0, gravity, 0);

        deltaq *= smoothingRadius;

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
        particleNeighbors = new NativeArray<int>(particleNumber * particleNumber, Allocator.Persistent);
        neighborCounter = new NativeArray<int>(particleNumber, Allocator.Persistent);

        gridManager = new GridManager(new Vector3[particleNumber], particleNumber, smoothingRadius);
    }

    public void Update()
    {
        float startTime = Time.realtimeSinceStartup;
        //add external forces
        addExternalForces();
        Debug.Log("The time passed for external force calculation " + (Time.realtimeSinceStartup - startTime) * 1000 + "ms");

        //Find neighbors
        startTime = Time.realtimeSinceStartup;
        if(FindNeighborsWithGrid)
        {
            findNeighborsOptimized();
        }
        else
        {
            findNeighbors();
        }
        Debug.Log("The time passed for finding neighboors " + (Time.realtimeSinceStartup - startTime) * 1000 + "ms");

        //calculate densities
        startTime = Time.realtimeSinceStartup;
        //calculateDensities();
        for (int i = 0; i < particleNumber; i++)
        {
            densities[i] = calculateDensity(i);
        }
        Debug.Log("The time passed for density calculation " + (Time.realtimeSinceStartup - startTime) * 1000 + "ms");

        startTime = Time.realtimeSinceStartup;
        for (int solveIteration = 0; solveIteration < solverIteration; solveIteration ++ )
        {
            //calculate lamda values
            calculateLamdas();

            //calculate deltaP
            calculateDeltaPs();

            //update predicted pos by respecting the constrain
            for (int i = 0; i < particleNumber; i++)
            {
                predictedPositions[i] += deltaPs[i]* deltaTime;
            }
        }
        Debug.Log("The time passed for constraint solver " + (Time.realtimeSinceStartup - startTime) * 1000 + "ms");

        //update position and velocity
        for (int i = 0; i < particleNumber; i++)
        {
            velocities[i] = (predictedPositions[i] - positions[i]) / deltaTime;
            positions[i] = predictedPositions[i];
            checkBoundaryConditions(i);
            particlesObjects[i].transform.position = positions[i];
        }

        clearNeighbors();

        debugNeighbors(0);
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

    private void addExternalForces()
    {
        //add gravity
        for (int i = 0; i < particleNumber; i++)
        {
            velocities[i] += getExternalForce() * deltaTime;
        }

        //add mouse interaction
        if(Input.GetMouseButton(0))
        {
            for(int i = 0; i < particleNumber; i++)
            {
                //get the mouse world position
                Vector3 mousePosition = Input.mousePosition;
                Vector3 worldPosition = Camera.main.ScreenToWorldPoint(new Vector3(mousePosition.x, mousePosition.y, 1));

                Vector3 interForce = InteractionForce(worldPosition, 15, 50, i) ;

                velocities[i] += new Vector3(interForce.x, interForce.y) * deltaTime;
            }
        }
        else if(Input.GetMouseButton(1))
        {
            for (int i = 0; i < particleNumber; i++)
            {
                //get the mouse world position
                Vector3 mousePosition = Input.mousePosition;
                Vector3 worldPosition = Camera.main.ScreenToWorldPoint(new Vector3(mousePosition.x, mousePosition.y, 1));

                Vector3 interForce = InteractionForce(worldPosition, 15, 50, i);

                velocities[i] -= new Vector3(interForce.x, interForce.y) * deltaTime;
            }
        }

        for (int i = 0; i < particleNumber; i++)
        {
            predictedPositions[i] = positions[i] + velocities[i] * deltaTime;
        }
    }
    private Vector3 InteractionForce(Vector3 inputPos, float radius, float strenght, int particleIndex)
    {
        //this function calculates the force applied to a particle to a specific position
        Vector3 interactionForce = Vector2.zero;
        Vector3 offset = inputPos - predictedPositions[particleIndex];
        float sqrDst = Vector2.Dot(offset, offset);

        if(sqrDst < radius * radius)
        {
            float dst = math.sqrt(sqrDst);
            Vector3 dirToInputPoint = dst <= float.Epsilon ? Vector2.zero : offset / dst;

            float centerT = 1 - dst / radius;

            interactionForce += (dirToInputPoint * strenght - velocities[particleIndex]) * centerT;
        }

        return interactionForce;
    }

    private void calculateDeltaPs()
    {

        // Create and schedule the job
        CalculateDeltaPs job = new CalculateDeltaPs
        {
            predictedPositions = predictedPositions,
            particleNeighbors = particleNeighbors, // assuming this is already a NativeArray
            neighborCounter = neighborCounter,
            smoothingRadius = smoothingRadius,
            particleNumber = particleNumber,
            deltaPs = deltaPs,
            lamdas = lamdas,
            activateSurfaceTension = activateSurfaceTension,
            deltaq = deltaq,
            k = k,
        };
        JobHandle handle = job.Schedule(particleNumber, 32); // Each job has equal number of particle to handle

        // Wait for the job to complete
        handle.Complete();

        /*for(int i = 0; i < particleNumber; i++)
        {
            deltaPs[i] = calculateDeltaP(i);
        }*/
    }

    private Vector3 calculateDeltaP(int i)
    {
        Vector3 deltaPi = Vector3.zero;
        float dist;
        Vector3 dir;
        int pj;

        for (int j = 0; j < neighborCounter[i]; j++)
        {
            pj = particleNeighbors[i * particleNumber + j];
            if (i == pj)
                continue;

            dist = (predictedPositions[i] - predictedPositions[pj]).magnitude;
            dir = dist == 0 ? getRandomDirection() : (predictedPositions[pj] - predictedPositions[i]) / dist;

            deltaPi += (lamdas[i] + lamdas[pj]) * gradientSpikyKernel(dist, smoothingRadius) * dir;
        }
        return deltaPi;
    }

    private void calculateLamdas()
    {
        CalculateLamdas job = new CalculateLamdas
        {
            predictedPositions = predictedPositions,
            particleNeighbors = particleNeighbors, // assuming this is already a NativeArray
            neighborCounter = neighborCounter,
            densities = densities,
            smoothingRadius = smoothingRadius,
            particleNumber = particleNumber,
            lamdas = lamdas,
            density = density,
            constraintForce = constraintForce,
        };
        JobHandle handle = job.Schedule(particleNumber, 32); // Each job has equal number of particle to handle

        // Wait for the job to complete
        handle.Complete();
    }

    private float calculateLamda(int i)
    {
        float denominator = 0;
        int k;

        for (int j = 0; j < neighborCounter[i]; j++)
        {
            k = particleNeighbors[i * particleNumber + j];
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
        int pj;

        if(k != i)
        {
            dist = (predictedPositions[i] - predictedPositions[k]).magnitude;
            partialDeriv = -gradientSpikyKernel(dist, smoothingRadius);
        }
        else
        {
            for (int j = 0; j < neighborCounter[i]; j++)
            {
                pj = particleNeighbors[i * particleNumber + j];
                dist = (predictedPositions[i] - predictedPositions[pj]).magnitude;
                partialDeriv += gradientSpikyKernel(dist, smoothingRadius);
            }
        }

        return partialDeriv / density;
    }

    private void calculateDensities()
    {
        // Create and schedule the job
        CalculateDensityJob job = new CalculateDensityJob
        {
            predictedPositions = predictedPositions,
            particleNeighbors = particleNeighbors, // assuming this is already a NativeArray
            neighborCounter = neighborCounter,
            particleMass = particleMass,
            smoothingRadius = smoothingRadius,
            densities = densities,
            particleNumber = particleNumber
        };
        JobHandle handle = job.Schedule(particleNumber, particleNumber/9); // Each job has equal number of particle to handle

        // Wait for the job to complete
        handle.Complete();
    }

    private float calculateDensity(int pi)
    {
        float density = 0;
        float dist;
        int pj;
        for(int j = 0; j < neighborCounter[pi]; j++)
        {
            pj = particleNeighbors[pi * particleNumber + j];
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
        if ((positions[pi].y - particleRadius) < verticalBoundaries.x && velocities[pi].y < 0)
        {
            velocities[pi] = new Vector3(velocities[pi].x, -velocities[pi].y * collisionDamping, velocities[pi].z);
            positions[pi] = new Vector3(positions[pi].x, verticalBoundaries.x, positions[pi].z);
        }
        if ((positions[pi].y + particleRadius) > verticalBoundaries.y && velocities[pi].y > 0)
        {
            velocities[pi] = new Vector3(velocities[pi].x, -velocities[pi].y * collisionDamping, velocities[pi].z);
            positions[pi] = new Vector3(positions[pi].x, verticalBoundaries.y, positions[pi].z);
        }
        //horizontal conditions
        if ((positions[pi].x - particleRadius) < horizontalBoundaries.x && velocities[pi].x < 0)
        {
            velocities[pi] = new Vector3(-velocities[pi].x * collisionDamping, velocities[pi].y, velocities[pi].z);
            positions[pi] = new Vector3(horizontalBoundaries.x, positions[pi].y, positions[pi].z);
        }
        if ((positions[pi].x + particleRadius) > horizontalBoundaries.y && velocities[pi].x > 0)
        {
            velocities[pi] = new Vector3(-velocities[pi].x * collisionDamping, velocities[pi].y, velocities[pi].z);
            positions[pi] = new Vector3(horizontalBoundaries.y, positions[pi].y, positions[pi].z);
        }
    }

    private void findNeighbors()
    {
        for(int i = 0; i < particleNumber; i++)
        {
            for(int j = 0; j < particleNumber; j++)
            {
                if ((predictedPositions[i] - predictedPositions[j]).magnitude < smoothingRadius)
                {
                    particleNeighbors[i*particleNumber + neighborCounter[i]] = j;
                    neighborCounter[i] += 1;
                }
            }
        }
    }

    private void findNeighborsOptimized()
    {
        gridManager.updatePositions(predictedPositions.ToArray());
        gridManager.findNeighbors(updateNeighbors);
    }
    private void updateNeighbors(int i, int j)
    {
        particleNeighbors[i * particleNumber + neighborCounter[i]] = j;
        neighborCounter[i] += 1;
    }

    private Vector3 getExternalForce()
    {
        return gravityVector;
    }

    private void clearNeighbors()
    {
        for (int i = 0; i < particleNumber; i++)
        {
            neighborCounter[i] = 0;
        }
    }

    private Vector3 getRandomDirection()
    {
        return new Vector3(UnityEngine.Random.Range(-1.0f, 1.0f), UnityEngine.Random.Range(-1.0f, 1.0f), 0).normalized;
    }

    private void debugNeighbors(int pi)
    {
        for(int j = 0; j < neighborCounter[pi]; j++)
        {
            int neighborIndex = particleNeighbors[pi * particleNumber + j];
            particlesObjects[neighborIndex].GetComponent<SpriteRenderer>().color = Color.red;
        }
    }

    void OnDestroy()
    {
        // Dispose NativeArrays to deallocate memory
        positions.Dispose();
        velocities.Dispose();
        lamdas.Dispose();
        deltaPs.Dispose();
        densities.Dispose();
        predictedPositions.Dispose();
        neighborCounter.Dispose();
        particleNeighbors.Dispose();
    }
}

struct CalculateDensityJob : IJobParallelFor
{
    [ReadOnly] public NativeArray<Vector3> predictedPositions; // NativeArray to store predicted positions
    [ReadOnly] public NativeArray<int> particleNeighbors;    
    [ReadOnly] public NativeArray<int> neighborCounter;  
    [ReadOnly] public float particleMass;
    [ReadOnly] public float smoothingRadius;
    public NativeArray<float> densities;            // NativeArray to store calculated densities
    [ReadOnly] public int particleNumber;

    public void Execute(int index)
    {
        float density = 0;
        float dist;

        // Loop through particle neighbors
        for (int j = 0; j < particleNeighbors.Length; j++)
        {
            int pj = particleNeighbors[j];
            dist = (predictedPositions[index] - predictedPositions[pj]).sqrMagnitude; // Use sqrMagnitude for squared distance
            density += particleMass * poly6Kernel(dist, smoothingRadius);
        }

        // Store the calculated density for this particle
        densities[index] = density;
    }

    // Your poly6Kernel function
    private float poly6Kernel(float r, float h)
    {
        if (r >= 0 && r <= h)
        {
            return (315 / (64 * math.PI * math.pow(h, 9))) * math.pow(h * h - r * r, 3);
        }
        return 0;
    }
}

struct CalculateDeltaPs : IJobParallelFor
{
    [ReadOnly] public NativeArray<Vector3> predictedPositions; // NativeArray to store predicted positions
    [ReadOnly] public NativeArray<int> particleNeighbors;
    [ReadOnly] public NativeArray<int> neighborCounter;
    [ReadOnly] public float smoothingRadius;
    [ReadOnly] public int particleNumber;
    public NativeArray<Vector3> deltaPs;
    [ReadOnly] public NativeArray<float> lamdas;
    [ReadOnly] public bool activateSurfaceTension;
    [ReadOnly] public float k;
    [ReadOnly] public float deltaq;

    public void Execute(int index)
    {
        Vector3 deltaPi = Vector3.zero;
        float dist;
        Vector3 dir;
        int pj;
        float Scorr = 0; //For Tensile Instability

        for (int j = 0; j < neighborCounter[index]; j++)
        {
            pj = particleNeighbors[index * particleNumber + j];
            if (index == pj)
                continue;

            dist = (predictedPositions[index] - predictedPositions[pj]).magnitude;
            dir = dist == 0 ? getRandomDirection() : (predictedPositions[pj] - predictedPositions[index]) / dist;

            if (activateSurfaceTension)
            {
                Scorr = -k * math.pow(poly6Kernel(dist, smoothingRadius) / poly6Kernel(deltaq, smoothingRadius), 4);
                //Debug.Log("Scorr is " + Scorr);
            }

            deltaPi += (lamdas[index] + lamdas[pj] + Scorr) * gradientSpikyKernel(dist, smoothingRadius) * dir;
        }
        deltaPs[index] = deltaPi;
    }

    private Vector3 getRandomDirection()
    {
        return new Vector3(0,1,0).normalized;
    }
    private float gradientSpikyKernel(float r, float h)
    {
        if (r >= 0 && r <= h)
        {
            return -(45 / (math.PI * math.pow(h, 6))) * math.pow(h - r, 2);
        }
        return 0;
    }

    private float poly6Kernel(float r, float h)
    {
        if (r >= 0 && r <= h)
        {
            return (315 / (64 * math.PI * math.pow(h, 9))) * math.pow(h * h - r * r, 3);
        }
        return 0;
    }
}

struct CalculateLamdas : IJobParallelFor
{
    [ReadOnly] public NativeArray<Vector3> predictedPositions; // NativeArray to store predicted positions
    [ReadOnly] public NativeArray<int> particleNeighbors;
    [ReadOnly] public NativeArray<int> neighborCounter;
    [ReadOnly] public NativeArray<float> densities;
    [ReadOnly] public float smoothingRadius;
    [ReadOnly] public int particleNumber;
    public NativeArray<float> lamdas;
    public float density;
    public float constraintForce;

    public void Execute(int index)
    {
        float denominator = 0;
        int k;

        for (int j = 0; j < neighborCounter[index]; j++)
        {
            k = particleNeighbors[index * particleNumber + j];
            denominator += math.pow(math.abs(constraintGradient(index, k)), 2);
        }

        lamdas[index] = -constraintFunction(index) / (denominator + constraintForce);
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
        int pj;

        if (k != i)
        {
            dist = (predictedPositions[i] - predictedPositions[k]).magnitude;
            partialDeriv = -gradientSpikyKernel(dist, smoothingRadius);
        }
        else
        {
            for (int j = 0; j < neighborCounter[i]; j++)
            {
                pj = particleNeighbors[i * particleNumber + j];
                dist = (predictedPositions[i] - predictedPositions[pj]).magnitude;
                partialDeriv += gradientSpikyKernel(dist, smoothingRadius);
            }
        }

        return partialDeriv / density;
    }

    private float gradientSpikyKernel(float r, float h)
    {
        if (r >= 0 && r <= h)
        {
            return -(45 / (math.PI * math.pow(h, 6))) * math.pow(h - r, 2);
        }
        return 0;
    }
}
