using UnityEngine;
using System;
using System.Runtime.InteropServices;

public class Bevy : MonoBehaviour {
	// Parameters visible in Unity interface
	public int xSteps = 64;
    public int ySteps = 64;
	public int count = 100;
		
    void Start () {
		// Create supershapes		
		for(int s = 0; s < count; s++)
		{
            // make game object and move it by a random Vector3
            GameObject newGameObject = new GameObject();
            newGameObject.transform.localScale *= 10;
            newGameObject.transform.Translate(150 * UnityEngine.Random.insideUnitSphere);
            
            // add the supershape script to the game object
            Supershape newShape = newGameObject.AddComponent<Supershape>();

            // give supershape random parameters in reasonable range
            newShape.a_a  = 1.0f;
            newShape.a_b  = 1.0f;
            newShape.a_m  = UnityEngine.Random.Range(1.0f, 8.0f);
            newShape.a_n1 = UnityEngine.Random.Range(1.0f, 5.0f);
            newShape.a_n2 = UnityEngine.Random.Range(1.0f, 5.0f);
            newShape.a_n3 = UnityEngine.Random.Range(1.0f, 5.0f);
            
            newShape.b_a  = 1.0f;
            newShape.b_b  = 1.0f;
            newShape.b_m  = UnityEngine.Random.Range(1.0f, 16.0f);
            newShape.b_n1 = UnityEngine.Random.Range(1.0f, 5.0f);
            newShape.b_n2 = UnityEngine.Random.Range(1.0f, 5.0f);
            newShape.b_n3 = UnityEngine.Random.Range(1.0f, 5.0f);

		}
		
        //    Debug.Log("Game object has no renderer or gui texture to assign the generated texture to!");
    }
    
    void OnDisable() {

    }
    
    void Update () {

    }
}
