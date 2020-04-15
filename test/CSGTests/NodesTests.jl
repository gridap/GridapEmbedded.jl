module NodesTests

using Test
using AbstractTrees
using GridapEmbedded.CSG

a = Leaf("a")
b = Leaf("b")
c = Leaf("c")

d = Node("d",a,b)
e = Node("e",c,d)
f = Node("f",e,b)
g = Node("g",f)

#print_tree(stdout,f)
#print_tree(stdout,g)

h = replace_data(objectid,f)
for (fi,hi) in zip(PreOrderDFS(f),PreOrderDFS(h))
  @test objectid(fi.data) == hi.data
end

h = replace_data(objectid,g)
for (gi,hi) in zip(PreOrderDFS(g),PreOrderDFS(h))
  @test objectid(gi.data) == hi.data
end

h = replace_data(d->d,objectid,f)
for (fi,hi) in zip(PreOrderDFS(f),PreOrderDFS(h))
  if isa(fi,Leaf)
    @test objectid(fi.data) == hi.data
  else
    @test fi.data === hi.data
  end
end

h = replace_data(d->d,objectid,g)
for (gi,hi) in zip(PreOrderDFS(g),PreOrderDFS(h))
  if isa(gi,Leaf)
    @test objectid(gi.data) == hi.data
  else
    @test gi.data === hi.data
  end
end

leaves = collect(Leaves(f))

leaves = collect(Leaves(g))

end # module
