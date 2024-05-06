using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class GridTest : MonoBehaviour
{
    //A test class to see if the Grid neighbor search works

    public GameObject particleObject;
    public float particleRadius = 0.5f;

    public int particleNumber = 200;
    public Vector2 cellNumber = Vector2.one;
    public Vector2 startPos = Vector2.zero;
    public float cellSize = 2;
    public float circleRadius = 2;

    private Vector3[] positions;

    private List<GameObject> particlesObjects = new List<GameObject>();

    private Vector2 horizontalBoundaries = new Vector2(0, 10);
    private Vector2 verticalBoundaries = new Vector2(0, 10);

    private GridManager gridManager;

    public void Start()
    {
        positions = new Vector3[particleNumber];

        //set the boundaries as camera boundaries
        Vector3 cameraCenter = Camera.main.transform.position;
        float halfHeight = Camera.main.orthographicSize;
        float halfWidth = Camera.main.aspect * halfHeight;
        verticalBoundaries = new Vector2(cameraCenter.y - halfHeight, cameraCenter.y + halfHeight);
        horizontalBoundaries = new Vector2(cameraCenter.x - halfWidth, cameraCenter.x + halfWidth);

        spawnRandomParticles();

        gridManager = new GridManager(positions, particleNumber, circleRadius);
    }

    public void Update()
    {
        reinitializeParticles();

        //get the mouse world position
        Vector3 mousePosition = Input.mousePosition;
        Vector3 worldPosition = Camera.main.ScreenToWorldPoint(new Vector3(mousePosition.x, mousePosition.y, 1));

        //find the neighbors and put their color to red
        gridManager.UpdateSpatialLookUp();
        gridManager.findNeighborPointsFromPosition(worldPosition, setParticleColorToRed);
    }

    void OnDrawGizmos()
    {
        drawGrid();

        Vector3 mousePosition = Input.mousePosition;
        Vector3 worldPosition = Camera.main.ScreenToWorldPoint(new Vector3(mousePosition.x, mousePosition.y, 1));
        Gizmos.DrawSphere(worldPosition, circleRadius);
    }

    private void drawGrid()
    {
        Vector3 verticalUp, center, horizontalRight;
        for (int i = 0; i < cellNumber.x; i++)
        {
            for (int j = 0; j < cellNumber.y; j++)
            {
                center = new Vector2(i * cellSize, j * cellSize);
                verticalUp = new Vector3(center.x, center.y + cellSize);
                horizontalRight = new Vector3(center.x + cellSize, center.y);

                Gizmos.DrawLine(center, verticalUp);
                Gizmos.DrawLine(center, horizontalRight);
            }
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
            positions[i] = randomSpawnPos;
        }
    }

    //this function puts all of the particle colors to blue
    public void reinitializeParticles()
    {
        for (int i = 0; i < particleNumber; i++)
        {
            particlesObjects[i].GetComponent<SpriteRenderer>().color = Color.blue;
        }
    }

    public void setParticleColorToRed(int i)
    {
        particlesObjects[i].GetComponent<SpriteRenderer>().color = Color.red;
    }
}
