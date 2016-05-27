using mesh
include("matrix.jl")
include("vector.jl")

#number of degrees of freedom in each direction
n1 = 20;
n2 = 20;

#Size of the rectangle
a = 1;
b = .5;

@printf "Size of the rectangle: %f %f 
Number of elements in each direction: %d %d \n" a b n1 n2 

m = mesh.Mesh(a,b,n1,n2);
mD = mesh.setBC(m);

#mesh.plotMeshWithBC(mD);
intorder = 1
K = globalMatrix(mD, intorder);

#pressure
h = 0.000005

#Right hand side zero
f(x) = zeros(size(x[:,1])[1],1)

B = globalVector(mD,f, h,intorder)

u = K\B

#plot displacements
mesh.plotDisplacements(mD,u)

