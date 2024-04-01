using System.Collections;
using System.Collections.Generic;
using Unity.Collections;
using Unity.Mathematics;
using Unity.VisualScripting;
using UnityEngine;

public class PBF_GPU : MonoBehaviour
{
    //This function handles the data to use it on a compute shader to execute the PBF simulation on the GPU
    public GameObject particleObject;

    [Header("Simulation parameters")]
    public int particleNumber;
    public Vector2 particleDistances = new Vector2(1, 1);
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
    //Parameters for viscocity
    public float viscocityCoef = 0.01f;
    public bool activateViscocity = false;

    [Header("References")]
    public ComputeShader compute;


    private List<GameObject> particlesObjects = new List<GameObject>();

    private Vector2 horizontalBoundaries = new Vector2(0, 10);
    private Vector2 verticalBoundaries = new Vector2(0, 10);
    private float2[] positions;
    private float2[] velocities;

    //Buffers
    ComputeBuffer positionsBuffer;
    ComputeBuffer predictedPositionsBuffer;
    ComputeBuffer velocitiesBuffer;
    ComputeBuffer densitiesBuffer;
    ComputeBuffer lambdasBuffer;
    ComputeBuffer deltaPsBuffer;

    ComputeBuffer spatialOffsetBuffer; //contains the sorted hash values corresponding the indice in the spacialIndice array
    ComputeBuffer spatialIndiciesBuffer; //Gives the start position of the particle in the same block

    enum SpawnType { UNIFORM, RANDOM };

    public void Start()
    {
        Debug.Log(compute.FindKernel("ExternalForces"));
        Debug.Log(compute.FindKernel("UpdatePositions"));

        positions = new float2[particleNumber];
        velocities = new float2[particleNumber];

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


        // Create buffers
        positionsBuffer = ComputeHelper.CreateStructuredBuffer<float2>(particleNumber);
        predictedPositionsBuffer = ComputeHelper.CreateStructuredBuffer<float2>(particleNumber);
        velocitiesBuffer = ComputeHelper.CreateStructuredBuffer<float2>(particleNumber);

        setBufferData();

        //Set the data to the buffer
        ComputeHelper.SetBuffer(compute, positionsBuffer, "Positions", 0, 1);
        ComputeHelper.SetBuffer(compute, predictedPositionsBuffer, "PredictedPositions", 0, 1);
        ComputeHelper.SetBuffer(compute, velocitiesBuffer, "Velocities", 0, 1);

        /*
        compute.SetBuffer(0, "Positions", positionsBuffer);
        compute.SetBuffer(0, "Velocities", velocitiesBuffer);
        compute.SetBuffer(0, "PredictedPositions", predictedPositionsBuffer);
        */

        //Set the constants to the buffer
        UpdateSetting(); 
    }

    public void Update()
    {
        //Debug.Log(string.Join(" ", positions));
        runSimulationStep();

        //get back the position
        positionsBuffer.GetData(positions);

        velocitiesBuffer.GetData(velocities);
        Debug.Log(string.Join(" ", positions));
        

    }

    public void LateUpdate()
    {
        //make the objects move
        for (int i = 0; i < particleNumber; i++)
        {
            particlesObjects[i].transform.position = new Vector3(positions[i].x, positions[i].y, 0);
        }
    }


    public void runSimulationStep()
    {
        ComputeHelper.Dispatch(compute, particleNumber, kernelIndex: 0);
        ComputeHelper.Dispatch(compute, particleNumber, kernelIndex: 1);
    }

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
            velocities[i] = float2.zero;
        }
    }

    public void UpdateSetting()
    {
        compute.SetFloat("gravity", gravity);
        compute.SetFloat("numParticle", particleNumber);
        compute.SetFloat("deltaTime", deltaTime);
    }

    public void setBufferData()
    {
        positionsBuffer.SetData(positions);
        velocitiesBuffer.SetData(velocities);
    }

    void OnDestroy()
    {
        ComputeHelper.Release(positionsBuffer, predictedPositionsBuffer, velocitiesBuffer, densitiesBuffer, spatialIndiciesBuffer, spatialOffsetBuffer);
    }
}
