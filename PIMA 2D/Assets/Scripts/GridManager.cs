using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Net;
using UnityEngine;
using UnityEngine.Experimental.GlobalIllumination;

public class GridManager : MonoBehaviour
{
    //This class proposes a Grid structure to structurize particle Positions. It contains some methode to draw and test the Grid. And also Useful
    //methods to analyze particles more efficiently

    public GameObject particleObject;
    public float particleRadius = 0.5f;

    public int particleNumber = 200;
    public Vector2 cellNumber = Vector2.one;
    public Vector2 startPos = Vector2.zero;
    public float cellSize = 2;
    public float circleRadius = 2;

    private List<(uint,uint)> spatialLookUp = new List<(uint, uint)>(); //first value = pos index second = hashval
    private int[] startIndices;
    private Vector2[] positions;

    private Vector2 horizontalBoundaries = new Vector2(0, 10);
    private Vector2 verticalBoundaries = new Vector2(0, 10);

    private List<(float, float)> cellOffsets = new List<(float, float)>();

    private List<GameObject> particlesObjects = new List<GameObject>();

    public void Start()
    {
        positions = new Vector2[particleNumber];
        startIndices = new int[particleNumber];

        //set the boundaries as camera boundaries
        Vector3 cameraCenter = Camera.main.transform.position;
        float halfHeight = Camera.main.orthographicSize;
        float halfWidth = Camera.main.aspect * halfHeight;
        verticalBoundaries = new Vector2(cameraCenter.y - halfHeight, cameraCenter.y + halfHeight);
        horizontalBoundaries = new Vector2(cameraCenter.x - halfWidth, cameraCenter.x + halfWidth);

        for (int i = -1; i < 2; i++)
        {
            for (int j = -1; j < 2; j++)
            {
                cellOffsets.Add((i * circleRadius, j * circleRadius));
            }
        }

        spawnRandomParticles();
    }

    public void Update()
    {
        UpdateSpatialLookUp(positions, circleRadius);

        //get the mouse world position
        Vector3 mousePosition = Input.mousePosition;
        Vector3 worldPosition = Camera.main.ScreenToWorldPoint(new Vector3(mousePosition.x, mousePosition.y, 1));

        findNeighborPoints(worldPosition);
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
        for(int i = 0; i < cellNumber.x; i ++)
        {
            for(int j = 0; j < cellNumber.y; j ++)
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

    //this is a tes function that find the neigbor positions and change their colors
    public void findNeighborPoints(Vector2 samplePoint)
    {
        (int centerX, int centerY) = PositionToCellCoord(samplePoint, circleRadius);

        //loop over all neighbor cells
        foreach((float offsetX, float offsetY) in cellOffsets)
        {
            //get the neighbor cells key
            (int cellIndexX, int cellIndexY) = PositionToCellCoord(new Vector2(samplePoint.x + offsetX, samplePoint.y + offsetY), circleRadius);
            Debug.Log((cellIndexX, cellIndexY));
            uint key = getKeyFromHash(HashCell(cellIndexX, cellIndexY), particleNumber);
            int cellStartIndex = startIndices[key];

            for (int i = cellStartIndex; i < spatialLookUp.Count; i ++)
            {
                //exit the loop if we are no longer looking at the correct cell
                if (spatialLookUp[i].Item2 != key) break;
             
                uint particleIndex = spatialLookUp[i].Item1;
                float dist = (positions[particleIndex] - samplePoint).magnitude;

                //check if the position is inside the radius
                
                particlesObjects[(int)particleIndex].GetComponent<SpriteRenderer>().color = Color.red;
            }
        }
    }

    //this function stocks the cell position for all points
    public void UpdateSpatialLookUp(Vector2[] points, float radius)
    {
        spatialLookUp.Clear();

        for(uint i = 0; i < points.Length; i++)
        {
            (int cellX, int cellY) = PositionToCellCoord(points[i], radius);
            uint cellKey = getKeyFromHash(HashCell(cellX, cellY), points.Length);
            spatialLookUp.Add((i, cellKey));
            startIndices[i] = 0;
        }

        //Sort by the cellKey
        spatialLookUp = spatialLookUp.OrderBy(pair => pair.Item2).ToList();

        //Calculate the start indices of each unique cell key in the spatial lookup
        for(int i = 0; i < points.Length; i ++)
        {
            uint key = spatialLookUp[i].Item2;
            uint keyPrev = i == 0 ? uint.MaxValue : spatialLookUp[i - 1].Item2;
            if(key != keyPrev)
            {
                startIndices[key] = i;
            }
        }

        /*
        Debug.Log("spatial lookup array: " + string.Join(", ", spatialLookUp));
        Debug.Log("startIndices array: " + string.Join(", ", startIndices));
        */
    }

    //Convert a position to the coordinate of the corresponding cell
    public (int x, int y) PositionToCellCoord(Vector2 point, float radius)
    {
        int cellX = (int)(point.x / radius);
        int cellY = (int)(point.y / radius);
        return (cellX, cellY);
    }

    //Convert a cell coordinate into a single number
    public uint HashCell(int cellX, int cellY)
    {
        //we choose 2 arbitrary prime numbers
        uint a = (uint)cellX * 262501;
        uint b = (uint)cellY * 457393;
        return a + b;
    }

    public uint getKeyFromHash(uint hash, int pointCount)
    {
        return hash % (uint)pointCount; 
    }

}
