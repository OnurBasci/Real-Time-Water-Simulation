using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Net;
using Unity.Collections;
using UnityEngine;
using UnityEngine.Experimental.GlobalIllumination;

public class GridManager
{
    //This class proposes a Grid structure to structurize particle Positions. This methode is usable from other classes

    public float particleRadius = 0.5f;

    public int particleNumber = 200;
    public float circleRadius;

    private List<(uint,uint)> spatialLookUp = new List<(uint, uint)>(); //first value = pos index second = hashval
    private int[] startIndices;
    private Vector3[] positions;

    private List<(float, float)> cellOffsets = new List<(float, float)>();

    public delegate void MyFunction(int particleIndex);
    public delegate void MyFunction2(int particleIndex, int neighborIndex);

    public GridManager(Vector3[] positions, int particleNumber, float circleRadius)
    {
        Debug.Log(circleRadius);
        this.positions = positions;
        this.particleNumber = particleNumber;
        this.circleRadius = circleRadius;

        startIndices = new int[particleNumber];


        for (int i = -1; i < 2; i++)
        {
            for (int j = -1; j < 2; j++)
            {
                cellOffsets.Add((i * circleRadius, j * circleRadius));
            }
        }
    }

    public void findNeighbors(MyFunction2 callBackFunction)
    {
        //This methods search adjacent neighbors of each position and apply the callbackfunctio if the neigbor is in the radius
        UpdateSpatialLookUp();
        
        for(int i = 0; i < particleNumber; i ++)
        {
            findNeighborPointsFromIndex(i, callBackFunction);
        }
    }

    //this is a tes function that find the neigbor positions and update neighbors. A delegetefunction is given to execute if the neighbors are found
    public void findNeighborPointsFromIndex(int particleIndex, MyFunction2 callBackFunction)
    {
        Vector3 pointPos = positions[particleIndex];
        //loop over all neighbor cells
        foreach ((float offsetX, float offsetY) in cellOffsets)
        {
            //get the neighbor cells key
            (int cellIndexX, int cellIndexY) = PositionToCellCoord(new Vector2(pointPos.x + offsetX, pointPos.y + offsetY), circleRadius);
            uint key = getKeyFromHash(HashCell(cellIndexX, cellIndexY), particleNumber);
            int cellStartIndex = startIndices[key];

            for (int i = cellStartIndex; i < spatialLookUp.Count; i++)
            {
                //exit the loop if we are no longer looking at the correct cell
                if (spatialLookUp[i].Item2 != key) break;

                uint neighborIndex = spatialLookUp[i].Item1;
                float dist = (positions[(int)neighborIndex] - pointPos).magnitude;
                //check if the position is inside the radius
                if (dist < circleRadius)
                {
                    //call the function to execute if the particle is a neighbor
                    callBackFunction(particleIndex, (int)neighborIndex);
                }
            }
        }
    }

    //this is a tes function that find the neigbor positions and change their colors. A delegetefunction is given to execute if the neighbors are found
    public void findNeighborPointsFromPosition(Vector3 samplePoint, MyFunction callBackFunction)
    {

        //loop over all neighbor cells
        foreach((float offsetX, float offsetY) in cellOffsets)
        {
            //get the neighbor cells key
            (int cellIndexX, int cellIndexY) = PositionToCellCoord(new Vector2(samplePoint.x + offsetX, samplePoint.y + offsetY), circleRadius);
            uint key = getKeyFromHash(HashCell(cellIndexX, cellIndexY), particleNumber);
            int cellStartIndex = startIndices[key];

            for (int i = cellStartIndex; i < spatialLookUp.Count; i ++)
            {
                //exit the loop if we are no longer looking at the correct cell
                if (spatialLookUp[i].Item2 != key) break;
             
                uint particleIndex = spatialLookUp[i].Item1;
                float dist = (positions[(int)particleIndex] - samplePoint).magnitude;
                Debug.Log(particleIndex);
                //check if the position is inside the radius
                if(dist < circleRadius)
                {
                    //call the function to execute if the particle is a neighbor

                }
                callBackFunction((int)particleIndex);
            }
        }
    }

    //this function stocks the cell position for all points
    public void UpdateSpatialLookUp()
    {
        spatialLookUp.Clear();

        for(uint i = 0; i < positions.Length; i++)
        {
            (int cellX, int cellY) = PositionToCellCoord(positions[i], circleRadius);
            uint cellKey = getKeyFromHash(HashCell(cellX, cellY), positions.Length);
            spatialLookUp.Add((i, cellKey));
            startIndices[i] = 0;
        }

        //Sort by the cellKey
        spatialLookUp = spatialLookUp.OrderBy(pair => pair.Item2).ToList();

        //Calculate the start indices of each unique cell key in the spatial lookup
        for(int i = 0; i < positions.Length; i ++)
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

    public void updatePositions(NativeArray<Vector3> pos)
    {
        for(int i = 0; i < particleNumber; i ++)
        {
            positions[i] = pos[i];
        }
    }
    public void updatePositions(Vector3[] pos)
    {
        positions = pos;
    }
}
