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

#print_tree(stdout,f)

h = replace_data(objectid,f)
for (fi,hi) in zip(PreOrderDFS(f),PreOrderDFS(h))
  @test objectid(fi.data) == hi.data
end

h = replace_data(d->d,objectid,f)
for (fi,hi) in zip(PreOrderDFS(f),PreOrderDFS(h))
  if isa(fi,Leaf)
    @test objectid(fi.data) == hi.data
  else
    @test fi.data === hi.data
  end
end

leaves = collect(Leaves(f))

end # module
