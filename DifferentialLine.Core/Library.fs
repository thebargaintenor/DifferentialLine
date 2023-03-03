namespace DifferentialLine.Core

open System.Collections.Generic

module Say =
    let hello name = printfn "Hello %s" name

type Node = { x: float; y: float }

type DifferentialLine =
    { nodes: List<Node>
      maxForce: float
      maxSpeed: float
      desiredSeparation: float
      cohesionRatio: float
      maxEdgeLength: float }

    member this.Add(n: Node) = this.nodes.Add(n)
    member this.Insert(n: Node) = this.nodes.Insert(n)
