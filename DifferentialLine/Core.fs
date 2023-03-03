namespace DifferentialLine.Core

open MathNet.Numerics.LinearAlgebra

module Core =
    type Node = { x: float; y: float }

    type DifferentialLine =
        { nodes: Node list
          maxForce: float
          maxSpeed: float
          desiredSeparation: float
          cohesionRatio: float
          maxEdgeLength: float }

    // this will have lots of opportunities for optimization
    let separationForces (dl: DifferentialLine) : Vector<double> list =
        let n = dl.nodes.Length
        let forces = Array.init n (fun i -> DenseVector.create 2 0.0)
        let nearNodes = Array.zeroCreate<int> n

        [DenseVector.create 2 0.0]




    let cohesionForces (dl: DifferentialLine) : Vector<double> list =
        dl.nodes |> List.map (fun n -> DenseVector.create 2 0.0)

    /// <summary>Calculates the changes in acceleration for each node
    /// in the given differential line path.
    ///
    /// For all nodes:
    /// dA/dt := Fseparation * scalingFactor + Fcohesion
    /// </summary>
    /// <returns>List of 2D acceleration deltas</returns>
    let differentiate (dl: DifferentialLine) : Vector<double> list =
        separationForces dl
        |> List.map (fun v -> v.Multiply dl.cohesionRatio)
        |> List.zip (cohesionForces dl)
        |> List.map (fun (a, b) -> a.Add(b))
