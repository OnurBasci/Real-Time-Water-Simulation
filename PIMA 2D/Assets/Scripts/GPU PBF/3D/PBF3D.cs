using System.Collections;
using System.Collections.Generic;
using Unity.Mathematics;
using UnityEngine;

public class PBF3D : MonoBehaviour
{
    //This function handles the data to use it on a compute shader to execute the PBF simulation on the GPU
    public GameObject particleObject;

    [Header("Simulation parameters")]
    public int particleNumber;
    public Vector3 boundaryLengths = new Vector3(5,10,15);
    public Vector3 particleDistances = new Vector3(1, 1, 1);
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
    //tension stability coefficients
    public float k = 0.1f; //a small positif coefficient
    public float deltaq = 0.1f;
    //Parameters for viscocity
    public float viscocityCoef = 0.01f;
    [Header("Mouse Interraction")]
    public float interactionStrength = 1;
    public float interactionRadius = 5;

    [Header("References")]
    public ComputeShader compute;


    private List<GameObject> particlesObjects = new List<GameObject>();

    private Vector3 horizontalBoundaries = new Vector2(0, 10);
    private Vector3 verticalBoundaries = new Vector2(0, 10);
    private Vector3 depthBoundaries = new Vector2(0, 10);
    private float3[] positions;
    private float3[] velocities;
    private uint[] spatialOffsets;
    private uint3[] spatialIndex;
    private uint[] particleNeighbors;
    private float[] densities;
    private float[] lamdas;
    private float3[] deltaPs;

    //Buffers
    ComputeBuffer positionsBuffer;
    ComputeBuffer predictedPositionsBuffer;
    ComputeBuffer velocitiesBuffer;
    ComputeBuffer densitiesBuffer;
    ComputeBuffer lambdasBuffer;
    ComputeBuffer deltaPsBuffer;
    ComputeBuffer particleNeighborBuffer;
    ComputeBuffer neighborCounterBuffer;
    GPUSort gpuSort;

    // Kernel IDs
    const int externalForcesKernel = 0;
    const int spatialHashKernel = 1;
    const int findNeighborKernel = 2;
    const int densityKernel = 3;
    const int lamdaKenrel = 4;
    const int deltaPKernel = 5;
    const int updatePredictedKernel = 6;
    const int updatePositionKernel = 7;

    ComputeBuffer spatialOffsetBuffer; //contains the sorted hash values corresponding the indice in the spacialIndice array
    ComputeBuffer spatialIndiciesBuffer; //Gives the start position of the particle in the same block

    enum SpawnType { UNIFORM, RANDOM };

    public void Start()
    {
        positions = new float3[particleNumber];
        velocities = new float3[particleNumber];
        spatialOffsets = new uint[particleNumber];
        spatialIndex = new uint3[particleNumber];
        particleNeighbors = new uint[particleNumber*particleNumber];
        densities = new float[particleNumber];
        lamdas = new float[particleNumber];
        deltaPs = new float3[particleNumber];

        deltaq *= smoothingRadius;

        //set the boundaries as camera boundaries
        horizontalBoundaries = new Vector2(-boundaryLengths.x/2, boundaryLengths.x/2);
        verticalBoundaries = new Vector2(-boundaryLengths.y / 2, boundaryLengths.y / 2);
        depthBoundaries = new Vector2(-boundaryLengths.z / 2, boundaryLengths.z / 2);

        /*if (spawnType == SpawnType.UNIFORM)
            spawnParticles();*/
        spawnRandomParticles();


        // Create buffers
        positionsBuffer = ComputeHelper.CreateStructuredBuffer<float3>(particleNumber);
        predictedPositionsBuffer = ComputeHelper.CreateStructuredBuffer<float3>(particleNumber);
        velocitiesBuffer = ComputeHelper.CreateStructuredBuffer<float3>(particleNumber);
        densitiesBuffer = ComputeHelper.CreateStructuredBuffer<float>(particleNumber);
        spatialIndiciesBuffer = ComputeHelper.CreateStructuredBuffer<uint3>(particleNumber);
        spatialOffsetBuffer = ComputeHelper.CreateStructuredBuffer<uint>(particleNumber);
        particleNeighborBuffer = ComputeHelper.CreateStructuredBuffer<uint>(particleNumber * particleNumber);
        neighborCounterBuffer = ComputeHelper.CreateStructuredBuffer<uint>(particleNumber);
        lambdasBuffer = ComputeHelper.CreateStructuredBuffer<float>(particleNumber);
        deltaPsBuffer = ComputeHelper.CreateStructuredBuffer<float3>(particleNumber);

        setBufferData();

        //Set the data to the buffer
        ComputeHelper.SetBuffer(compute, positionsBuffer, "Positions", externalForcesKernel, updatePositionKernel);
        ComputeHelper.SetBuffer(compute, predictedPositionsBuffer, "PredictedPositions", externalForcesKernel, spatialHashKernel,
            updatePositionKernel, densityKernel, findNeighborKernel, lamdaKenrel, deltaPKernel, updatePredictedKernel);
        ComputeHelper.SetBuffer(compute, velocitiesBuffer, "Velocities", externalForcesKernel, updatePositionKernel);
        ComputeHelper.SetBuffer(compute, densitiesBuffer, "Densities", densityKernel, lamdaKenrel);
        ComputeHelper.SetBuffer(compute, spatialIndiciesBuffer, "SpatialIndices", spatialHashKernel, findNeighborKernel, densityKernel);
        ComputeHelper.SetBuffer(compute, spatialOffsetBuffer, "SpatialOffsets", spatialHashKernel, findNeighborKernel, densityKernel);
        ComputeHelper.SetBuffer(compute, particleNeighborBuffer, "particleNeighbors", densityKernel, updatePositionKernel, findNeighborKernel, lamdaKenrel, deltaPKernel);
        ComputeHelper.SetBuffer(compute, neighborCounterBuffer, "neighborCounter", findNeighborKernel, densityKernel, updatePositionKernel, lamdaKenrel, deltaPKernel);
        ComputeHelper.SetBuffer(compute, lambdasBuffer, "lamdas", lamdaKenrel, deltaPKernel);
        ComputeHelper.SetBuffer(compute, deltaPsBuffer, "deltaPs", deltaPKernel, updatePredictedKernel);

        //Set the constants to the buffer
        UpdateSetting();

        gpuSort = new();
        gpuSort.SetBuffers(spatialIndiciesBuffer, spatialOffsetBuffer);
    }

    public void Update()
    {
        //Check mouse interraction
        checkMouseInterraction();

        //Debug.Log(particleNumber);
        runSimulationStep();

        //get back the position
        positionsBuffer.GetData(positions);

        deltaPsBuffer.GetData(deltaPs);
        //Debug.Log("offsets " + string.Join(" ", deltaPs));

        //make the objects move
        //float startTime = Time.realtimeSinceStartup;
        for (int i = 0; i < particleNumber; i++)
        {
            particlesObjects[i].transform.position = new Vector3(positions[i].x, positions[i].y, positions[i].z);
        }
        //Debug.Log("The time passed for moving particles " + (Time.realtimeSinceStartup - startTime) * 1000 + "ms");
    }



    public void runSimulationStep()
    {
        ComputeHelper.Dispatch(compute, particleNumber, kernelIndex: externalForcesKernel);
        ComputeHelper.Dispatch(compute, particleNumber, kernelIndex: spatialHashKernel);
        gpuSort.SortAndCalculateOffsets();
        ComputeHelper.Dispatch(compute, particleNumber, kernelIndex: findNeighborKernel);
        for (int simulationIter = 0; simulationIter < solverIteration; simulationIter++)
        {
            ComputeHelper.Dispatch(compute, particleNumber, kernelIndex: densityKernel);
            ComputeHelper.Dispatch(compute, particleNumber, kernelIndex: lamdaKenrel);
            ComputeHelper.Dispatch(compute, particleNumber, kernelIndex: deltaPKernel);
            ComputeHelper.Dispatch(compute, particleNumber, kernelIndex: updatePredictedKernel);
        }
        ComputeHelper.Dispatch(compute, particleNumber, kernelIndex: updatePositionKernel);
    }

    /*
    private void spawnParticles()
    {

        Vector3 cameraCenter = Camera.main.transform.position;
        float2 spawnPos;

        int numColumn = (int)Mathf.Sqrt(particleNumber);
        int lastLineParticles = particleNumber % numColumn;

        float2 offset = new float2(numColumn / 2 * particleDistances.x, numColumn / 2 * particleDistances.y);

        for (int i = 0; i < particleNumber; i++)
        {
            int row = i / numColumn;
            int column = i % numColumn;

            spawnPos = new float2(cameraCenter.x + column * particleDistances.x, cameraCenter.y + row * particleDistances.y) - offset;
            GameObject p1 = Instantiate(particleObject, new Vector3(spawnPos.x, spawnPos.y, 0), Quaternion.identity);
            p1.transform.localScale = new Vector3(particleRadius * 2, particleRadius * 2, particleRadius * 2);
            particlesObjects.Add(p1);
            positions[i] = spawnPos;
            velocities[i] = float2.zero;
        }

    }
    */

    public void spawnRandomParticles()
    {
        Vector3 randomSpawnPos;
        for (int i = 0; i < particleNumber; i++)
        {

            randomSpawnPos = new Vector3(UnityEngine.Random.Range(horizontalBoundaries.x, horizontalBoundaries.y),
                UnityEngine.Random.Range(verticalBoundaries.x, verticalBoundaries.y), UnityEngine.Random.Range(depthBoundaries.x, depthBoundaries.y));
            GameObject p1 = Instantiate(particleObject, randomSpawnPos, Quaternion.identity);
            p1.transform.localScale = new Vector3(particleRadius * 2, particleRadius * 2, particleRadius * 2);
            particlesObjects.Add(p1);

            //intialize particle data
            positions[i] = randomSpawnPos;
            velocities[i] = float3.zero;
        }
    }

    public void UpdateSetting()
    {
        compute.SetFloat("gravity", gravity);
        compute.SetInt("numParticle", particleNumber);
        compute.SetFloat("deltaTime", deltaTime);
        compute.SetFloat("collisionDamping", collisionDamping);
        compute.SetFloat("verticalBoundariesLower", verticalBoundaries.x);
        compute.SetFloat("verticalBoundariesUpper", verticalBoundaries.y);
        compute.SetFloat("horizontalBoundariesLeft", horizontalBoundaries.x);
        compute.SetFloat("horizontalBoundariesRight", horizontalBoundaries.y);
        compute.SetFloat("depthBoundariesFront", depthBoundaries.x);
        compute.SetFloat("depthBoundariesBack", depthBoundaries.y);
        compute.SetFloat("particleRadius", particleRadius);
        compute.SetFloat("smoothingRadius", smoothingRadius);
        compute.SetFloat("targetDensity", density);
        compute.SetFloat("surfaceTensionParam1", k);
        compute.SetFloat("surfaceTensionParam2", deltaq);
        compute.SetFloat("viscocityCoef", viscocityCoef);
        compute.SetVector("interactionInputPoint", new Vector2(0, 0));
        compute.SetFloat("interactionInputStrength", 0);
        compute.SetFloat("interactionInputRadius", 0);
    }

    public void checkMouseInterraction()
    {
        // Mouse interaction settings:
        bool isPullInteraction = Input.GetMouseButton(0);
        bool isPushInteraction = Input.GetMouseButton(1);
        float currInteractStrength;
        if (isPushInteraction || isPullInteraction)
        {
            currInteractStrength = isPushInteraction ? -interactionStrength : interactionStrength;

            Vector2 mousePos = Camera.main.ScreenToWorldPoint(Input.mousePosition);

            compute.SetVector("interactionInputPoint", mousePos);
            compute.SetFloat("interactionInputStrength", currInteractStrength);
            compute.SetFloat("interactionInputRadius", interactionRadius);
        }
        else
        {
            compute.SetVector("interactionInputPoint", new Vector2(0, 0));
            compute.SetFloat("interactionInputStrength", 0);
            compute.SetFloat("interactionInputRadius", 0);
        }
    }

    public void setBufferData()
    {
        positionsBuffer.SetData(positions);
        velocitiesBuffer.SetData(velocities);
    }


    void OnDestroy()
    {
        ComputeHelper.Release(positionsBuffer, predictedPositionsBuffer, velocitiesBuffer, densitiesBuffer,
            spatialIndiciesBuffer, spatialOffsetBuffer, neighborCounterBuffer, particleNeighborBuffer, lambdasBuffer, deltaPsBuffer);
    }

    private void clearNeighbors()
    {
        ComputeHelper.Release(particleNeighborBuffer, neighborCounterBuffer);

        particleNeighborBuffer = ComputeHelper.CreateStructuredBuffer<int>(particleNumber * particleNumber);
        neighborCounterBuffer = ComputeHelper.CreateStructuredBuffer<int>(particleNumber);

        ComputeHelper.SetBuffer(compute, particleNeighborBuffer, "particleNeighbors", findNeighborKernel);
        ComputeHelper.SetBuffer(compute, neighborCounterBuffer, "neighborCounter", findNeighborKernel);
    }

    void OnDrawGizmosSelected()
    {
        // Draw a cube to show the borders
        Gizmos.color = new Color(0.2f, 0.2f, 0.2f, 0.5f);
        Gizmos.DrawCube(Vector3.zero, new Vector3(boundaryLengths.x, boundaryLengths.y, boundaryLengths.z));
    }
}
