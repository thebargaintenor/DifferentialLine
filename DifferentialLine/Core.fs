namespace DifferentialLine.Core

open MathNet.Numerics
open MathNet.Numerics.LinearAlgebra

module VecMath =
    let scaleToMagnitude (magnitude: float) (vector: Vector<float>): Vector<float> =
        vector.Multiply(magnitude / vector.L2Norm())

    let add (a: Vector<'a>) (b: Vector<'a>): Vector<'a> = b.Add(a)

    let subtract (a: Vector<'a>) (b: Vector<'a>): Vector<'a> = b.Subtract(a)

    let multiplyScalar (a: 'a) (b: Vector<'a>): Vector<'a> = b.Multiply(a)

    let divideScalar (a: 'a) (b: Vector<'a>): Vector<'a> = b.Divide(a)

    let limit (max: float) (vector: Vector<float>): Vector<float> =
        if vector.L2Norm() > max then
            scaleToMagnitude max vector
        else
            vector
    
    /// Euclidean space unit vectors are calculated with 2-norm
    let Normalize (vector: Vector<float>) = vector.Normalize(2)

module Core =
    type Node = 
        { position: Vector<float>; 
          velocity: Vector<float>; 
          fmax: float;
          vmax: float; }

    type DifferentialLine =
        { nodes: Node array
          maxForce: float
          maxSpeed: float
          desiredSeparation: float
          cohesionRatio: float
          maxEdgeLength: float }

    let origin = DenseVector.create 2 0.0

    let separationForce (cutoffDistance: float) (a: Node) (b: Node): Vector<float> =
        let dist = Distance.Euclidean(a.position, b.position)
        if dist > 0 && dist <= cutoffDistance then
            a.position.Subtract(b.position)
            |> VecMath.Normalize
            |> VecMath.divideScalar dist
        else
            DenseVector.create 2 0.0

    // this will have lots of opportunities for optimization
    let separationForces (dl: DifferentialLine) : Vector<float> array =
        let n = dl.nodes.Length
        let forces = Array.init n (fun i -> origin)
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
                        VecMath.scaleToMagnitude dl.maxSpeed forces.[i]
                        |> VecMath.subtract dl.nodes.[i].velocity
                        |> VecMath.limit dl.maxForce
        forces

    let attract (node: Node) (target: Vector<float>): Vector<float> =
        target
        |> VecMath.subtract node.position
        |> VecMath.scaleToMagnitude node.vmax
        |> VecMath.subtract node.velocity
        |> VecMath.limit node.fmax


    let cohesionForces (dl: DifferentialLine) : Vector<float> array =
        let nmax = dl.nodes.Length - 1
        
        [| 0 .. nmax |]
        |> Array.map (fun i ->
            match i with
            | 0 -> dl.nodes.[nmax].position.Add(dl.nodes.[1].position)
            | n when n = nmax -> dl.nodes.[n-1].position.Add(dl.nodes.[0].position)
            | _ -> dl.nodes.[i-1].position.Add(dl.nodes.[i+1].position)
            )
        |> Array.map (VecMath.divideScalar 2.0)
        |> Array.map2 attract dl.nodes

    /// <summary>Calculates the changes in acceleration for each node
    /// in the given differential line path.
    ///
    /// For all nodes:
    /// dA/dt := Fseparation * scalingFactor + Fcohesion
    /// </summary>
    /// <returns>List of 2D acceleration deltas</returns>
    let differentiate (dl: DifferentialLine) : Vector<float> array =
        separationForces dl
        |> Array.map (VecMath.multiplyScalar dl.cohesionRatio)
        |> Array.map2 VecMath.add (cohesionForces dl)
