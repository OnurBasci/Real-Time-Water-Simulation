using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.Rendering;

public class TextureManipulator : MonoBehaviour
{
    public GameObject squareObject;

    public float thresh = 0;

    private Texture2D bgTexture;

    private GameObject backgroundObject;

    public int width = 256;
    public int height = 256;
    public SPHSimulation sphSim;

    private SPHSimulation.Particle[] particles;
    private float[] colorValues;

    float widthOfBG, heightOfBG;

    Renderer bgRenderer;

    private void Start()
    {

        initializeBackgroundObject();

        bgTexture = new Texture2D(width, height);

        bgTexture.wrapMode = TextureWrapMode.Clamp;
        bgTexture.filterMode = FilterMode.Bilinear;
    }

    private void OnGUI()
    {
        GUI.depth = -1;
        GUI.DrawTexture(new Rect(0f, 0f, Screen.width, Screen.height), bgTexture);
    }

    void Update()
    {
        particles = sphSim.particles;
        colorValues = sphSim.densities;
        UpdateTextureFromArray();
    }

    private void initializeBackgroundObject()
    {
        Vector3 cameraCenter = Camera.main.transform.position;
        float halfHeight = Camera.main.orthographicSize;
        float halfWidth = Camera.main.aspect * halfHeight;

        width = (int) (height*Camera.main.aspect);

        widthOfBG = halfWidth * 2;
        heightOfBG = halfHeight * 2;

        //backgroundObject = Instantiate(squareObject, cameraCenter, Quaternion.identity);

        //backgroundObject.transform.localScale = new Vector3(widthOfBG, heightOfBG, 0);
    }

    void UpdateTextureFromArray()
    {

        /*Debug.Log(bgTexture.width);

        for (int y = 0; y < bgTexture.height; y++)
        {
            for (int x = 0; x < bgTexture.width; x++)
            {
                Debug.Log("(" + x + ", " + y +")");
                
                bgTexture.SetPixel(x, y, Color.red);
                
            }
        }

        bgTexture.Apply();*/


        //bgRenderer.material.SetTexture("_MainTex", bgTexture);
        //bgRenderer.material.mainTexture = bgTexture;



        Color pixelColor;
        for(int i = 0; i < colorValues.Length; i ++)
        {
            
            int x = (int)(width * particles[i].position.x / widthOfBG);
            int y = (int)(height * particles[i].position.y / heightOfBG);

            if (colorValues[i] > thresh)
            {
                bgTexture.SetPixel(x, y, Color.red);
            }
            else
            {
                bgTexture.SetPixel(x, y, Color.blue);
            }
        }

        bgTexture.Apply();
    }
}
