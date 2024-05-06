using System.Collections;
using System.Collections.Generic;
using System.Linq;
using UnityEngine;

public class Obstacle : MonoBehaviour
{
    public enum ObstacleType
    {
        CIRCLE,
        RECTANGLE
    }
    public FluidShower fluidShower;
    public ObstacleType obstacleType;

    private HashSet<GameObject> blocks = new HashSet<GameObject>();
    private HashSet<GameObject> oldblocks = new HashSet<GameObject>();
    private HashSet<GameObject> freedObjects;

    private Rigidbody2D rb;

    private void Awake()
    {
        rb = GetComponent<Rigidbody2D>();
    }


    public void Update()
    {
        getCollidedBlocks();
        fluidShower.fluid.setBlocksAsObsticles(blocks, rb.velocity);

        freedObjects = new HashSet<GameObject>(oldblocks);
        freedObjects.Except(blocks);

        Debug.Log(freedObjects.Count);

        fluidShower.fluid.setBlocksAsNonObsticles(freedObjects);

        oldblocks = blocks;
        blocks.Clear();

    }

    private void getCollidedBlocks()
    {
        Collider2D[] colliders = null;
        if(obstacleType == ObstacleType.CIRCLE)
        {
            float radius = GetComponent<CircleCollider2D>().radius;
            colliders = Physics2D.OverlapCircleAll(transform.position, radius);
        }
        else if(obstacleType == ObstacleType.RECTANGLE)
        {
            Bounds bounds = GetComponent<BoxCollider2D>().bounds;
            colliders = Physics2D.OverlapAreaAll(new Vector2(bounds.min.x, bounds.min.y), new Vector2(bounds.max.x, bounds.max.y));
        }
        foreach (Collider2D col in colliders)
        {
            blocks.Add(col.gameObject);
        }
    }
}
