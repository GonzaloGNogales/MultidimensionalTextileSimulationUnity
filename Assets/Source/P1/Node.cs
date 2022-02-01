﻿using System.Collections;
using System.Collections.Generic;
using UnityEngine;

using VectorXD = MathNet.Numerics.LinearAlgebra.Vector<double>;
using MatrixXD = MathNet.Numerics.LinearAlgebra.Matrix<double>;
using DenseVectorXD = MathNet.Numerics.LinearAlgebra.Double.DenseVector;
using DenseMatrixXD = MathNet.Numerics.LinearAlgebra.Double.DenseMatrix;

public class Node : MonoBehaviour {

    #region InEditorVariables

    public float Mass;
    public bool Fixed;

    #endregion

    public int index;

    public Vector3 Pos;
    public Vector3 Vel;

    private PhysicsManager Manager;

	// Update is called once per frame
	void Update () {
        transform.position = Pos;
	}

    // Use this for initialization
    public void Initialize(int ind, PhysicsManager m)
    {
        index = ind;
        Manager = m;

        Pos = transform.position;
    }

    public void GetPosition(VectorXD pos)
    {
        pos[index] = Pos.x;
        pos[index + 1] = Pos.y;
        pos[index + 2] = Pos.z;
    }

    public void SetPosition(VectorXD pos)
    {
        Pos = new Vector3((float)pos[index], (float)pos[index + 1], (float)pos[index + 2]);
    }

    public void GetVelocity(VectorXD vel)
    {
        vel[index] = Vel.x;
        vel[index + 1] = Vel.y;
        vel[index + 2] = Vel.z;
    }

    public void SetVelocity(VectorXD vel)
    {
        Vel = new Vector3((float)vel[index], (float)vel[index + 1], (float)vel[index + 2]);
    }

    public void GetForce(VectorXD force)
    {
        // Add gravity and damping forces related with actual node vel
        
        // Spring Force
        Vector3 Force = Mass * Manager.Gravity;
        // Damping Force
        Force += - 0.4f * Mass * Vel;
        
        force[index] += Force.x;
        force[index + 1] += Force.y;
        force[index + 2] += Force.z;
    }

    // Get Force Jacobian
    public void GetForceJacobian(MatrixXD dFdx, MatrixXD dFdv)
    {
        // The derivative of (0, m*g.y, 0) is => dFdx = (0, 0, 0) as there is no x in the force (only 2nd Newton's Law)
        // But we have to manage dFdv for simulating damping force
        // Fill dFdv (D) matrix
        // Row 0 dFadxa
        // dFdv[index, index] += -0.4f * Mass;
        // dFdv[index, index + 1] += -0.4f * Mass;
        // dFdv[index, index + 2] += -0.4f * Mass;
        // // Row 1 dFadxa
        dFdv[index + 1, index] += -0.4f * Mass;
        dFdv[index + 1, index + 1] += -0.4f * Mass;
        dFdv[index + 1, index + 2] += -0.4f * Mass;
        // // Row 2 dFadxa
        // dFdv[index + 2, index] += -0.4f * Mass;
        // dFdv[index + 2, index + 1] += -0.4f * Mass;
        // dFdv[index + 2, index + 2] += -0.4f * Mass;
    }

    public void GetMass(MatrixXD mass)
    {
        mass[index, index] = Mass;
        mass[index+1, index+1] = Mass;
        mass[index+2, index+2] = Mass;
    }

    public void GetMassInverse(MatrixXD massInv)
    {
        massInv[index, index] = 1.0 / Mass;
        massInv[index + 1, index + 1] = 1.0 / Mass;
        massInv[index + 2, index + 2] = 1.0 / Mass;
    }

    public void FixVector(VectorXD v)
    {
        if(Fixed)
        {
            v[index] = 0.0;
            v[index+1] = 0.0;
            v[index+2] = 0.0;
        }
    }

    public void FixMatrix(MatrixXD M)
    {
        if (Fixed)
        {
            for (int i = 0; i < M.RowCount; i++)
            {
                M[index, i] = 0.0;
                M[index+1, i] = 0.0;
                M[index+2, i] = 0.0;
                M[i, index] = 0.0;
                M[i, index+1] = 0.0;
                M[i, index+2] = 0.0;
            }
            M[index, index] = 1.0;
            M[index+1, index+1] = 1.0;
            M[index+2, index+2] = 1.0;
        }
    }

}
