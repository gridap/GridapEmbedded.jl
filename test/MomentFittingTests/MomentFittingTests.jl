module MomentFittingTest

  using Test
  using Gridap
  using GridapEmbedded
  using GridapEmbedded.MomentFitting

  function run_test(domain,partition,geom,degree)
    bgmodel = CartesianDiscreteModel(domain,partition)
    cutgeo = cut(bgmodel,geom)
    Ωᵃ = Triangulation(cutgeo,ACTIVE_IN,geom)
    dΩᵃ = Measure(MomentFittingQuad(Ωᵃ,cutgeo,degree))
    Ωᶜ = Triangulation(cutgeo)
    dΩᶜ = Measure(Ωᶜ,num_dims(bgmodel)*degree,degree)
    @test sum(∫(1)dΩᵃ) - sum(∫(1)dΩᶜ) + 1 ≈ 1
    if num_dims(bgmodel) == 2
      uᵃ = CellField(x->(x[1]+x[2])^degree,Ωᵃ)
      uᶜ = CellField(x->(x[1]+x[2])^degree,Ωᶜ)
    else
      uᵃ = CellField(x->(x[1]+x[2]+x[3])^degree,Ωᵃ)
      uᶜ = CellField(x->(x[1]+x[2]+x[3])^degree,Ωᶜ)
    end
    @test sum(∫(uᵃ)dΩᵃ) - sum(∫(uᶜ)dΩᶜ) + 1 ≈ 1
  end

  # CASE 2D 1
  n = 6
  partition = (n,n)
  domain = (0,1,0,1)
  geom = disk(0.42,x0=Point(0.5,0.5))
  # println("CASE 2D 1")
  for deg in 1:10 # 19
    # println(deg)
    run_test(domain,partition,geom,deg)
  end

  # CASE 2D 2
  n = 6
  partition = (n,n)
  domain = (0,1,0,1)
  geom = square(L=0.53,x0=Point(0.5,0.5))
  # println("CASE 2D 2")
  for deg in 1:10 # 19
    # println(deg)
    run_test(domain,partition,geom,deg)
  end

  # CASE 2D 3
  n = 6
  m = 7
  partition = (n,m)
  domain = (0,1,0,1)
  geom = disk(0.42,x0=Point(0.5,0.5))
  # println("CASE 2D 3")
  for deg in 1:10 # 19
    # println(deg)
    run_test(domain,partition,geom,deg)
  end

  # CASE 3D 1
  n = 6
  partition = (n,n,n)
  domain = (0,1,0,1,0,1)
  geom = sphere(0.42,x0=Point(0.5,0.5,0.5))
  # println("CASE 3D 1")
  for deg in 1:5 # 5
    # println(deg)
    run_test(domain,partition,geom,deg)
  end

  # CASE 3D 2
  n = 6
  partition = (n,n,n)
  domain = (0,1,0,1,0,1)
  geom = cube(L=0.53,x0=Point(0.5,0.5,0.5))
  # println("CASE 3D 2")
  for deg in 1:5 # 7
    # println(deg)
    run_test(domain,partition,geom,deg)
  end

  # CASE 3D 3
  n = 2
  partition = (n,n,n)
  domain = (0,1,0,1,0,1)
  geom = plane(x0=Point(0.3,0.5,0.5),v=VectorValue(-1.0,0.0,0.0))
  # println("CASE 3D 3")
  for deg in 1:5 # 11
    # println(deg)
    run_test(domain,partition,geom,deg)
  end

end # module

