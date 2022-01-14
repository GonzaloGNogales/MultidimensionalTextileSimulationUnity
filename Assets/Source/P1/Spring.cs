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
        // Direction of the Forces
        MatrixXD u = new DenseMatrixXD(1, 3);
        Vector3 dir = nodeA.Pos - nodeB.Pos;
        dir.Normalize();
        u[0, 0] = dir[0];
        u[0, 1] = dir[1];
        u[0, 2] = dir[2];
        
        // Identity matrix
        MatrixXD I = DenseMatrixXD.CreateIdentity(3);
        
        MatrixXD dFadxa = - Stiffness * (Length - Length0) / Length * (I - u.Transpose() * u) - Stiffness * u.Transpose() * u;
        MatrixXD dFadxb = - dFadxa;
        MatrixXD dFbdxa = - dFadxa;
        MatrixXD dFbdxb = dFadxa;
        
        // Fill dFdx (K) matrix
        // Row 0 dFadxa
        dFdx[nodeA.index, nodeA.index] += dFadxa[0, 0];
        dFdx[nodeA.index, nodeA.index + 1] += dFadxa[0, 1];
        dFdx[nodeA.index, nodeA.index + 2] += dFadxa[0, 2];
        // Row 1 dFadxa
        dFdx[nodeA.index + 1, nodeA.index] += dFadxa[1, 0];
        dFdx[nodeA.index + 1, nodeA.index + 1] += dFadxa[1, 1];
        dFdx[nodeA.index + 1, nodeA.index + 2] += dFadxa[1, 2];
        // Row 2 dFadxa
        dFdx[nodeA.index + 2, nodeA.index] += dFadxa[2, 0];
        dFdx[nodeA.index + 2, nodeA.index + 1] += dFadxa[2, 1];
        dFdx[nodeA.index + 2, nodeA.index + 2] += dFadxa[2, 2];
        
        // Row 0 dFadxb
        dFdx[nodeA.index, nodeB.index] += dFadxb[0, 0];
        dFdx[nodeA.index, nodeB.index + 1] += dFadxb[0, 1];
        dFdx[nodeA.index, nodeB.index + 2] += dFadxb[0, 2];
        // Row 1 dFadxb
        dFdx[nodeA.index + 1, nodeB.index] += dFadxb[1, 0];
        dFdx[nodeA.index + 1, nodeB.index + 1] += dFadxb[1, 1];
        dFdx[nodeA.index + 1, nodeB.index + 2] += dFadxb[1, 2];
        // Row 2 dFadxb
        dFdx[nodeA.index + 2, nodeB.index] += dFadxb[2, 0];
        dFdx[nodeA.index + 2, nodeB.index + 1] += dFadxb[2, 1];
        dFdx[nodeA.index + 2, nodeB.index + 2] += dFadxb[2, 2];
        
        // Row 0 dFbdxa
        dFdx[nodeB.index, nodeA.index] += dFbdxa[0, 0];
        dFdx[nodeB.index, nodeA.index + 1] += dFbdxa[0, 1];
        dFdx[nodeB.index, nodeA.index + 2] += dFbdxa[0, 2];
        // Row 1 dFbdxa
        dFdx[nodeB.index + 1, nodeA.index] += dFbdxa[1, 0];
        dFdx[nodeB.index + 1, nodeA.index + 1] += dFbdxa[1, 1];
        dFdx[nodeB.index + 1, nodeA.index + 2] += dFbdxa[1, 2];
        // Row 2 dFbdxa
        dFdx[nodeB.index + 2, nodeA.index] += dFbdxa[2, 0];
        dFdx[nodeB.index + 2, nodeA.index + 1] += dFbdxa[2, 1];
        dFdx[nodeB.index + 2, nodeA.index + 2] += dFbdxa[2, 2];
        
        // Row 0 dFbdxb
        dFdx[nodeB.index, nodeB.index] += dFbdxb[0, 0];
        dFdx[nodeB.index, nodeB.index + 1] += dFbdxb[0, 1];
        dFdx[nodeB.index, nodeB.index + 2] += dFbdxb[0, 2];
        // Row 1 dFbdxb
        dFdx[nodeB.index + 1, nodeB.index] += dFbdxb[1, 0];
        dFdx[nodeB.index + 1, nodeB.index + 1] += dFbdxb[1, 1];
        dFdx[nodeB.index + 1, nodeB.index + 2] += dFbdxb[1, 2];
        // Row 2 dFbdxb
        dFdx[nodeB.index + 2, nodeB.index] += dFbdxb[2, 0];
        dFdx[nodeB.index + 2, nodeB.index + 1] += dFbdxb[2, 1];
        dFdx[nodeB.index + 2, nodeB.index + 2] += dFbdxb[2, 2];
    }

}
