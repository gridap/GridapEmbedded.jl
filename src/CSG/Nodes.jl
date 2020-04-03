
struct Node{Td,Ta,Tb}
  data::Td
  leftchild::Ta
  rightchild::Tb
end

Base.eltype(::Type{<:Node{Td}}) where Td = Td

const Leaf = Node{Td,Nothing,Nothing} where Td

function Leaf(data)
  Node(data,nothing,nothing)
end

function children(a::Node)
  (a.leftchild,a.rightchild)
end

function children(a::Leaf)
  ()
end

function Base.show(io::IO,a::Node)
  print(io,"Node(")
  show(io,a.data)
  print(io,",")
  show(io,a.leftchild)
  print(io,",")
  show(io,a.rightchild)
  print(io,")")
end

function Base.show(io::IO,a::Leaf)
  print(io,"Leaf(")
  show(io,a.data)
  print(io,")")
end

function replace_data(filter,a::Node)
  data = filter(a.data)
  c1 = replace_data(filter,a.leftchild)
  c2 = replace_data(filter,a.rightchild)
  Node(data,c1,c2)
end

function replace_data(filter,a::Leaf)
  data = filter(a.data)
  Leaf(data)
end

function replace_data(filter_node,filter_leaf,a::Node)
  data = filter_node(a.data)
  c1 = replace_data(filter_node,filter_leaf,a.leftchild)
  c2 = replace_data(filter_node,filter_leaf,a.rightchild)
  Node(data,c1,c2)
end

function replace_data(filter_node,filter_leaf,a::Leaf)
  data = filter_leaf(a.data)
  Leaf(data)
end

