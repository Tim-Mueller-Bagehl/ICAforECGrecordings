using Pkg
Pkg.activate(@__DIR__)
Pkg.resolve()
Pkg.instantiate()

using ICAforECGrecordings
Pkg.test("ICAforECGrecordings")