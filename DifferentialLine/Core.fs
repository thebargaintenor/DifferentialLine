namespace DifferentialLine.Core

open MathNet.Numerics
open MathNet.Numerics.LinearAlgebra

module VecMath =
    let ScaleToMagnitude (magnitude: float) (vector: Vector<float>): Vector<float> =
        vector.Multiply(magnitude / vector.L2Norm())

    let Add (a: Vector<'a>) (b: Vector<'a>): Vector<'a> = b.Add(a)

    let Subtract (a: Vector<'a>) (b: Vector<'a>): Vector<'a> = b.Subtract(a)

    let MultiplyScalar (a: 'a) (b: Vector<'a>): Vector<'a> = b.Multiply(a)

    let DivideScalar (a: 'a) (b: Vector<'a>): Vector<'a> = b.Divide(a)

    let Limit (max: float) (vector: Vector<float>): Vector<float> =
        if vector.L2Norm() > max then
            ScaleToMagnitude max vector
        else
            vector
    
    /// Euclidean space unit vectors are calculated with 2-norm
    let Normalize (vector: Vector<float>) = vector.Normalize(2)

module Core =
    type Node = 
        { position: Vector<float>; 
          velocity: Vector<float>; }

    type DifferentialLine =
        { nodes: Node array
          maxForce: float
          maxSpeed: float
          desiredSeparation: float
          cohesionRatio: float
          maxEdgeLength: float }

    let separationForce (cutoffDistance: float) (a: Node) (b: Node): Vector<float> =
        let dist = Distance.Euclidean(a.position, b.position)
        if dist > 0 && dist <= cutoffDistance then
            a.position.Subtract(b.position)
            |> VecMath.Normalize
            |> VecMath.DivideScalar dist
        else
            DenseVector.create 2 0.0

    // this will have lots of opportunities for optimization
    let separationForces (dl: DifferentialLine) : Vector<float> array =
        let n = dl.nodes.Length
        let forces = Array.init n (fun i -> DenseVector.create 2 0.0)
        let nearNodes = Array.zeroCreate<int> n

        for i = 0 to n-1 do
            for j = i+1 to n-1 do
                let repulsion = separationForce dl.desiredSeparation dl.nodes.[i] dl.nodes.[j]
                if repulsion.L2Norm() > 0 then
                    forces.[i] <- forces.[i].Add(repulsion)
                    forces.[j] <- forces.[j].Add(repulsion)
                    nearNodes.[i] <- nearNodes.[i] + 1
                    nearNodes.[j] <- nearNodes.[j] + 1

                if nearNodes.[i] > 0 then
                    forces.[i] <- forces.[i].Divide(nearNodes.[i])

                if forces.[i].L2Norm() > 0 then
                    forces.[i] <- 
                        VecMath.ScaleToMagnitude dl.maxSpeed forces.[i]
                        |> VecMath.Subtract dl.nodes.[i].velocity
                        |> VecMath.Limit dl.maxForce
        forces

    let cohesionForces (dl: DifferentialLine) : Vector<float> array =
        dl.nodes |> Array.map (fun n -> DenseVector.create 2 0.0)

    /// <summary>Calculates the changes in acceleration for each node
    /// in the given differential line path.
    ///
    /// For all nodes:
    /// dA/dt := Fseparation * scalingFactor + Fcohesion
    /// </summary>
    /// <returns>List of 2D acceleration deltas</returns>
    let differentiate (dl: DifferentialLine) : Vector<float> array =
        separationForces dl
        |> Array.map (fun v -> v.Multiply dl.cohesionRatio)
        |> Array.zip (cohesionForces dl)
        |> Array.map (fun (a, b) -> a.Add(b))
