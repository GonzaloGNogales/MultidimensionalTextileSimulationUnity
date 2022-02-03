using UnityEngine;
using System.Collections.Generic;
using VectorXD = MathNet.Numerics.LinearAlgebra.Vector<double>;
using MatrixXD = MathNet.Numerics.LinearAlgebra.Matrix<double>;
using DenseVectorXD = MathNet.Numerics.LinearAlgebra.Double.DenseVector;
using DenseMatrixXD = MathNet.Numerics.LinearAlgebra.Double.DenseMatrix;

/// <summary>
/// Basic physics manager capable of simulating a given ISimulable
/// implementation using diverse integration methods: explicit,
/// implicit, Verlet and semi-implicit.
/// </summary>
public class PhysicsManager : MonoBehaviour 
{
	/// <summary>
	/// Default constructor. Zero all. 
	/// </summary>
	public PhysicsManager()
	{
		Paused = true;
		TimeStep = 0.01f;
		Gravity = new Vector3 (0.0f, -9.81f, 0.0f);
		IntegrationMethod = Integration.Explicit;
        Substeps = 1;
    }

	/// <summary>
	/// Integration method.
	/// </summary>
	public enum Integration
	{
		Explicit = 0,
		Symplectic = 1,
        Midpoint = 2,
        Verlet = 3,
        Implicit = 4,
	};

	#region InEditorVariables

	public bool Paused;
	public float TimeStep;
    public Vector3 Gravity;
    public List<GameObject> SimObjects;
    public Integration IntegrationMethod;
    public int Substeps;
    
    #endregion

    #region OtherVariables
    private List<ISimulable> m_objs;
    private int m_numDoFs;
    
    // Verlet integrarion method variables
    private bool firstVerletStep;
    private VectorXD xOld;
    #endregion

    #region MonoBehaviour

    public void Start()
    {
        //Parse the simulable objects and initialize their state indices
        m_numDoFs = 0;
        m_objs = new List<ISimulable>(SimObjects.Capacity);

        foreach (GameObject obj in SimObjects)
        {
            ISimulable simobj = obj.GetComponent<ISimulable>();
            if (simobj != null)
            {
                m_objs.Add(simobj);

                // Initialize simulable object
                simobj.Initialize(m_numDoFs, this);

                // Retrieve pos and vel size
                m_numDoFs += simobj.GetNumDoFs();
            }
        }
        
        // Initialize Verlet integration variables
        firstVerletStep = true;
        xOld = new DenseVectorXD(m_numDoFs);
    }

    public void Update()
	{
		if (Input.GetKeyUp (KeyCode.P))
			Paused = !Paused;
    }

    public void FixedUpdate()
    {
        if (Paused)
            return; // Not simulating

        // Select integration method and add sub-stepping for speeding up the simulation
        for (int i = 0; i < Substeps; ++i)
        {
            switch (IntegrationMethod)
            {
                case Integration.Explicit: stepExplicit(); break;
                case Integration.Symplectic: stepSymplectic(); break;
                case Integration.Midpoint: stepMidpoint(); break;
                case Integration.Verlet: stepVerlet(); break;
                case Integration.Implicit: stepImplicit(); break;
                default:
                    throw new System.Exception("[ERROR] Should never happen!");
            }
        }
    }

    #endregion

    /// <summary>
    /// Performs a simulation step using Explicit integration.
    /// </summary>
    private void stepExplicit()
	{
        VectorXD x = new DenseVectorXD(m_numDoFs);
        VectorXD v = new DenseVectorXD(m_numDoFs);
        VectorXD f = new DenseVectorXD(m_numDoFs);
        MatrixXD Minv = new DenseMatrixXD(m_numDoFs);

        foreach (ISimulable obj in m_objs)
        {
            obj.GetPosition(x);
            obj.GetVelocity(v);
            obj.GetForce(f);
            obj.GetMassInverse(Minv);
        }

        foreach (ISimulable obj in m_objs)
        {
            obj.FixVector(f);
            obj.FixMatrix(Minv);
        }

        x += TimeStep * v;
        v += TimeStep * (Minv * f);

        foreach (ISimulable obj in m_objs)
        {
            obj.SetPosition(x);
            obj.SetVelocity(v);
        }
    }

    /// <summary>
    /// Performs a simulation step using Symplectic integration.
    /// </summary>
    private void stepSymplectic()
	{
        VectorXD x = new DenseVectorXD(m_numDoFs);
        VectorXD v = new DenseVectorXD(m_numDoFs);
        VectorXD f = new DenseVectorXD(m_numDoFs);
        MatrixXD Minv = new DenseMatrixXD(m_numDoFs);

        foreach (ISimulable obj in m_objs)
        {
            obj.GetPosition(x);
            obj.GetVelocity(v);
            obj.GetForce(f);
            obj.GetMassInverse(Minv);
        }

        foreach (ISimulable obj in m_objs)
        {
            obj.FixVector(f);
            obj.FixMatrix(Minv);
        }

        v += TimeStep * (Minv * f);
        x += TimeStep * v;

        foreach (ISimulable obj in m_objs)
        {
            obj.SetPosition(x);
            obj.SetVelocity(v);
        }
    }

    /// <summary>
    /// Performs a simulation step using Midpoint integration.
    /// </summary>
    private void stepMidpoint()
    {
        VectorXD x = new DenseVectorXD(m_numDoFs);
        VectorXD v = new DenseVectorXD(m_numDoFs);
        VectorXD f = new DenseVectorXD(m_numDoFs);
        MatrixXD Minv = new DenseMatrixXD(m_numDoFs);

        // Initialize and evaluate forces at t0
        foreach (ISimulable obj in m_objs)
        {
            obj.GetPosition(x);
            obj.GetVelocity(v);
            obj.GetForce(f);
            obj.GetMassInverse(Minv);
        }
        
        // Save initial position and velocity
        VectorXD x0 = x;
        VectorXD v0 = v;

        foreach (ISimulable obj in m_objs)
        {
            obj.FixVector(f);
            obj.FixMatrix(Minv);
        }

        // Midpoint Pos and Vel integration (with implicit position integration)
        v += 0.5f * TimeStep * (Minv * f);
        x += 0.5f * TimeStep * v;
        
        // After each integration it is needed to update spring state
        foreach (ISimulable obj in m_objs)
        {
            // This is called for updating spring states and
            // for computing forces again with correct mid positions and velocities
            obj.SetPosition(x);  
            obj.SetVelocity(v);
        }
        
        // Evaluate forces at h/2
        foreach (ISimulable obj in m_objs)
        { 
            // Perform another force computation
            obj.GetForce(f);
        }

        // Fix forces again because they may have changed
        foreach (ISimulable obj in m_objs)
        {
            obj.FixVector(f);
        }

        // Full step Pos and Vel integration (implicit Pos)
        v = v0 + TimeStep * (Minv * f);
        x = x0 +  TimeStep * v;
        
        // Finally set final positions and velocities with the midpoint integration resulting values of x and v
        foreach (ISimulable obj in m_objs)
        {
            obj.SetPosition(x);  // This is called for updating spring states and updating the positions
            obj.SetVelocity(v);
        }
    }

    /// <summary>
    /// Performs a simulation step using Verlet integration.
    /// </summary>
    private void stepVerlet()
    {
        VectorXD x = new DenseVectorXD(m_numDoFs);
        VectorXD v = new DenseVectorXD(m_numDoFs);
        VectorXD f = new DenseVectorXD(m_numDoFs);

        foreach (ISimulable obj in m_objs)
        {
            obj.GetPosition(x);
            obj.GetVelocity(v);
            obj.GetForce(f);
        }

        foreach (ISimulable obj in m_objs)
        {
            obj.FixVector(f);
        }

        if (firstVerletStep)  // Initialization supposing acceleration is constant
        {
            xOld = x;
            x += TimeStep * v + 0.5 * TimeStep * TimeStep * f;
            firstVerletStep = false;
        }
        else
        {
            VectorXD aux = x;
            x = 2 * x - xOld + TimeStep * TimeStep * f;
            v = (x - xOld) / (2 * TimeStep);
            xOld = aux;
        }

        foreach (ISimulable obj in m_objs)
        {
            obj.SetPosition(x);
            obj.SetVelocity(v);
        }
    }

    /// <summary>
    /// Performs a simulation step using Implicit integration.
    /// </summary>
    private void stepImplicit()
    {
        VectorXD x = new DenseVectorXD(m_numDoFs);
        VectorXD v = new DenseVectorXD(m_numDoFs);
        VectorXD f = new DenseVectorXD(m_numDoFs);
        VectorXD b = new DenseVectorXD(m_numDoFs);

        MatrixXD A = new DenseMatrixXD(m_numDoFs);
        MatrixXD M = new DenseMatrixXD(m_numDoFs);
        MatrixXD dFdx = new DenseMatrixXD(m_numDoFs);  // Elastic force derivative
        MatrixXD dFdv = new DenseMatrixXD(m_numDoFs);  // Damping force derivative
        
        foreach (ISimulable obj in m_objs)
        {
            obj.GetPosition(x);
            obj.GetVelocity(v);
            obj.GetForce(f);
            obj.GetMass(M);
            obj.GetForceJacobian(dFdx, dFdv);
        }
        
        /*
        // Check dFdx values outside of the foor loop for having dFdx complete
        // Check dFdx with finite differences: dFdx approx => (F(x0 + d) - F(x0)) / d
        float d = 0.001f;  // Approximate d checking with more than one value 0.01f and 0.001f for example
        VectorXD xSave = x;
        VectorXD fSave = f;
        //obj.GetPosition(x + d);
        //obj.GetForce(f);
        VectorXD approxdFdx = (f - fSave) / d;
        // Use approx_dFdx to check if derivatives were correctly computed
            
        // After checking restore state
        x = xSave;
        f = fSave;
        /**/

        foreach (ISimulable obj in m_objs)
        {
            obj.FixVector(f);
            obj.FixMatrix(M);
            obj.FixMatrix(dFdx);
            obj.FixMatrix(dFdv);
        }

        // The velocity implicit integration is computed solving a linear system
        A = M - TimeStep * dFdv - TimeStep * TimeStep * dFdx;
        b = (M - TimeStep * dFdv) * v + TimeStep * f;
        v = A.Solve(b);  // Solving A * v(t+h) = b => v(t+h)
        // The position integration in Multidimensional is the same as it was already implicit
        x += TimeStep * v;      

        foreach (ISimulable obj in m_objs)
        {
            obj.SetPosition(x);
            obj.SetVelocity(v);
        }
    }

}
