using System;
using System.Collections.Generic;
using Unity.Mathematics;
using UnityEngine;

public class SPHSimulation : MonoBehaviour
{
    public GameObject particleObject;

    public float gravity = -9.81f;
    public float density = 1;

    public int particleNumber = 100;
    public float particleRadius = 0.5f;
    public float particleHorizontalDistance = 2;
    public float particleVerticalDistance = 2;

    public float smoothingRadius = 5; //the radius for the smoothing function

    public float K_temp = 1; //coefficient depending on the temperature
    public float viscosityCoeff = 1;

    public float collisionDamping = 1;
    public Vector2 horizontalBoundaries = new Vector2(0, 10);
    public Vector2 verticalBoundaries = new Vector2(0, 10);

    [SerializeField] private SpawnType spawnType;

    private List<GameObject> particlesObjects = new List<GameObject>();
    [NonSerialized] public Particle[] particles;

    [NonSerialized] public float[] densities;
    private Vector3[] predictedPositions;

    enum SpawnType { UNIFORM, RANDOM };

    public class Particle
    {
        public Vector3 position;
        public Vector3 velocity;
        public float mass;
        public float radius;

        public Particle(Vector3 pos, Vector3 vel, float mass, float radius)
        {
            this.position = pos;
            this.velocity = vel;
            this.mass = mass;
            this.radius = radius;
        }
    }

    private void Start()
    {
        densities = new float[particleNumber];
        predictedPositions = new Vector3[particleNumber];

        //set the boundaries as camera boundaries
        Vector3 cameraCenter = Camera.main.transform.position;
        float halfHeight = Camera.main.orthographicSize;
        float halfWidth = Camera.main.aspect * halfHeight;
        verticalBoundaries = new Vector2(cameraCenter.y - halfHeight, cameraCenter.y + halfHeight);
        horizontalBoundaries = new Vector2(cameraCenter.x - halfWidth, cameraCenter.x + halfWidth);


        /*GameObject p1 = Instantiate(particleObject, cameraCenter, Quaternion.identity);
        particlesObjects.Add(p1);
        particles.Add(new Particle(new Vector3(5, 5, 5), new Vector3(0, 0, 0), 1, 0.5f));*/
        particles = new Particle[particleNumber];

        if (spawnType == SpawnType.UNIFORM)
            spawnParticles();
        else if (spawnType == SpawnType.RANDOM)
            spawnRandomParticles();
    }


    private void spawnParticles()
    {

        Vector3 cameraCenter = Camera.main.transform.position;
        Vector3 spawnPos;

        int numColumn = (int)Mathf.Sqrt(particleNumber);
        int lastLineParticles = particleNumber % numColumn;

        Vector3 offset = new Vector3(numColumn / 2 * particleHorizontalDistance, numColumn / 2 * particleVerticalDistance);

        for (int i = 0; i < particleNumber; i++)
        {
            int row = i / numColumn;
            int column = i % numColumn;

            spawnPos = new Vector3(cameraCenter.x + column * particleHorizontalDistance, cameraCenter.y + row * particleVerticalDistance, 0) - offset;
            GameObject p1 = Instantiate(particleObject, spawnPos, Quaternion.identity);
            p1.transform.localScale = new Vector3(particleRadius * 2, particleRadius * 2, particleRadius * 2);
            particlesObjects.Add(p1);
            particles[i] = new Particle(spawnPos, new Vector3(0, 0, 0), 1, particleRadius);
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
            particles[i] = new Particle(randomSpawnPos, new Vector3(0, 0, 0), 1, particleRadius);
        }
    }


    private void Update()
    {
        for(int i = 0; i <particleNumber; i++)
        {
            particles[i].velocity += addExternalForces() * Time.deltaTime;
        }

        //predict positions
        for (int i = 0; i < particleNumber; i++)
        {
            predictedPositions[i] = particles[i].position + particles[i].velocity * Time.deltaTime;
        }
        //calculate densities
        for (int i = 0; i < particleNumber; i++)
        {
            densities[i] = calculateDensity(predictedPositions[i]);
        }

        moveParticles();
    }

    private void moveParticles()
    {
        Particle p_i;
        for (int i = 0; i < particles.Length; i++)
        {
            //calculate the forces to apply
            p_i = particles[i];
            Vector3 appliedForces = calculateForce(i);

            //update the particles position
            Vector3 accelartion = appliedForces / densities[i];
            p_i.velocity += accelartion * Time.deltaTime;
            checkBoundaryConditions(p_i);
            p_i.position += p_i.velocity * Time.deltaTime;

            particlesObjects[i].transform.position = p_i.position;
        }
    }

    private void checkBoundaryConditions(Particle p)
    {
        //vertical conditions
        if ((p.position.y - p.radius) < verticalBoundaries.x && p.velocity.y < 0)
        {
            p.velocity.y = -p.velocity.y * collisionDamping;
        }
        if ((p.position.y + p.radius) > verticalBoundaries.y && p.velocity.y > 0)
        {
            p.velocity.y = -p.velocity.y * collisionDamping;
        }
        //horizontal conditions
        if ((p.position.x - p.radius) < horizontalBoundaries.x && p.velocity.x < 0)
        {
            p.velocity.x = -p.velocity.x * collisionDamping;
        }
        if ((p.position.x + p.radius) > horizontalBoundaries.y && p.velocity.x > 0)
        {
            p.velocity.x = -p.velocity.x * collisionDamping;
        }
    }

    private Vector3 calculateForce(int particleIndex)
    {
        //the force is calculated as the sum of the forces coming from pressure, viscocity and external forces
        Vector3 appliedForce = Vector3.zero;

        appliedForce += addPressurseForce(particleIndex);
        appliedForce += addViscocityForce(particleIndex);
        //appliedForce += addExternalForces();

        return appliedForce;
    }


    private Vector3 addPressurseForce(int particleIndex)
    {
        Vector3 f_pressure = Vector3.zero;
        Vector3 dir;
        float press_i = K_temp * (densities[particleIndex] - density);
        float press_j, dist;

        Particle pi = particles[particleIndex];
        Particle pj;

        for (int j = 0; j < particleNumber; j++)
        {
            if (j == particleIndex)
                continue;

            pj = particles[j];

            dist = (predictedPositions[particleIndex] - predictedPositions[j]).magnitude;

            dir = dist == 0 ? getRandomDirection() : (predictedPositions[j] - predictedPositions[particleIndex]) / dist;
            press_j = K_temp * (densities[j] - density);
            f_pressure += (pj.mass * (press_i + press_j) / (2 * densities[j])) * gradientSpikyKernel(dist, smoothingRadius) * dir;
        }
        return -f_pressure;
    }

    private Vector3 addViscocityForce(int particleIndex)
    {
        /*Vector3 f_viscocity = Vector3.zero;

        Particle pi = particles[particleIndex];
        Particle pj;

        Vector3 dir;
        float dist;

        for(int j = 0; j < particleNumber; j++)
        {
            if (j == particleIndex) continue;

            pj = particles[j];

            dist = (predictedPositions[particleIndex] - predictedPositions[j]).magnitude;

            dir = dist == 0 ? getRandomDirection() : (predictedPositions[j] - predictedPositions[particleIndex]) / dist;

            f_viscocity += pj.mass * ((pj.velocity.magnitude - pi.velocity.magnitude) / densities[j]) * laplaceSpikyKernel(dist, smoothingRadius) *dir;
        }
        */
        //return viscosityCoeff * f_viscocity;
        return Vector3.zero;
    }

    private Vector3 addExternalForces()
    {
        return new Vector3(0, density * gravity, 0);
    }

    private float calculateDensity(Vector3 r)
    {
        float density = 0;
        float dist;
        Particle pj;

        for(int j = 0; j < particles.Length; j ++)
        {
            pj = particles[j];
            dist = (r - predictedPositions[j]).magnitude;  //TO DO dist is supposed to be squared
            density += pj.mass * poly6Kernel(dist, smoothingRadius);
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

    private float laplaceSpikyKernel(float r, float h)
    {
        if(r >= 0 && r <= h)
        {
            return (90 / (math.PI * math.pow(h, 6))) * (h - r);
        }
        return 0;
    }

    private Vector3 getRandomDirection()
    {
        return new Vector3(UnityEngine.Random.Range(-1.0f, 1.0f), UnityEngine.Random.Range(-1.0f, 1.0f), 0).normalized;
    }

    void setPressureValues()
    {

    }
}
