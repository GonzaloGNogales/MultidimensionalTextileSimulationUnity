using System.Collections;
using System.Collections.Generic;
using UnityEngine;

using VectorXD = MathNet.Numerics.LinearAlgebra.Vector<double>;
using MatrixXD = MathNet.Numerics.LinearAlgebra.Matrix<double>;
using DenseVectorXD = MathNet.Numerics.LinearAlgebra.Double.DenseVector;
using DenseMatrixXD = MathNet.Numerics.LinearAlgebra.Double.DenseMatrix;

public class Spring : MonoBehaviour {

    #region InEditorVariables

    public float Stiffness;
    public Node nodeA;
    public Node nodeB;

    #endregion

    public float Length0;
    public float Length;

    private PhysicsManager Manager;

    // Update is called once per frame
    void Update () {
        Vector3 yaxis = new Vector3(0.0f, 1.0f, 0.0f);
        Vector3 dir = nodeA.Pos - nodeB.Pos;
        dir.Normalize();

        transform.position = 0.5f * (nodeA.Pos + nodeB.Pos);
        //The default length of a cylinder in Unity is 2.0
        transform.localScale = new Vector3(transform.localScale.x, Length / 2.0f, transform.localScale.z);
        transform.rotation = Quaternion.FromToRotation(yaxis, dir);
	}

    // Use this for initialization
    public void Initialize(PhysicsManager m)
    {
        Manager = m;

        UpdateState();
        Length0 = Length;
    }

    // Update spring state
    public void UpdateState()
    {
        Length = (nodeA.Pos - nodeB.Pos).magnitude;
    }

    // Get Force
    public void GetForce(VectorXD force)
    {
        // Add Hooke's law and damping forces related with actual nodes vel
        // Direction of the Forces
        Vector3 u = nodeA.Pos - nodeB.Pos;
        u.Normalize();
        // Spring Force
        Vector3 Force = - Stiffness * (Length - Length0) * u;
        // Damping Force
        Force += - 0.01f * Stiffness * Vector3.Dot(u, nodeA.Vel - nodeB.Vel) * u;
        
        // Node A
        force[nodeA.index] += Force.x;
        force[nodeA.index + 1] += Force.y;
        force[nodeA.index + 2] += Force.z;
        
        // Node B
        force[nodeB.index] -= Force.x;
        force[nodeB.index + 1] -= Force.y;
        force[nodeB.index + 2] -= Force.z;
    }

    // Get Force Jacobian
    public void GetForceJacobian(MatrixXD dFdx)
    {
        // TO BE COMPLETED //
    }

}
