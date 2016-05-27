module mesh

export Mesh, plotMesh, setDirichlet

Pkg.add("PyPlot")
using PyPlot

type Mesh
  width
  height
  coords
  elems::Array{Int,2}
  #Constructor
  function Mesh(width,height,n1,n2)
    @assert (width > 0 && height > 0 && n1 > 0 && n2 > 0) "Mesh arguments must be positive"
    coords = createCoords(width,height,n1,n2)
    elems = createElems(coords, n1,n2)
    new(width,height,coords,elems)
  end
end

function createCoords(a,b,n1,n2)
  #Fill coords array
  coords = zeros(2,(n1+1)*(n2+1));
  cnt = 1
  for i = 1:n2+1
    for j = 1:n1+1
      coords[1:2,cnt] = [a/n1*(j-1),b/n2*(i-1)]
      cnt = cnt+1
    end
  end
  return coords
end

#fill elems array
function createElems(coords, n1,n2)
  elems = zeros(4,n1*n2)
  cnt = 1
  for i = 1:n1
    for j = 1:n2
      elems[1:4,cnt] = [(j-1)*(n1+1)+i, (j-1)*(n1+1)+i+1, j*(n1+1)+i+1, j*(n1+1)+i]; 
      cnt  = cnt +1
    end
  end
  return elems
end

function plotMesh(m::Mesh)
  #Plot coordinates
  for i = 1:size(m.coords)[2]
    plot(m.coords[1,i],m.coords[2,i],marker="x",color ="k")
  end
  #plot elements
  for i = 1:size(m.elems)[2]
    plot(m.coords[1,m.elems[1:4,i]]',m.coords[2,m.elems[1:4,i]]',color = "k")
  end
end

dof = 2

type MeshWithBC
  mesh::Mesh
  ID::Array{Int,2}
  LM::Array{Int,2}
  IN::Array{Int,1}
  #Constructor
  function MeshWithBC(mesh,ID,LM,IN)
    new(mesh,ID,LM,IN)
  end
end

function plotMeshWithBC(mBC::MeshWithBC)
  m = mBC.mesh
  plotMesh(m)
  for i = 1:size(mBC.ID)[2]
    if mBC.ID[1,i] == 0
      plot(m.coords[1,i],m.coords[2,i],marker="o",color ="b");
    end
    if mBC.ID[1,i] in mBC.IN
       plot(m.coords[1,i],m.coords[2,i],marker="o",color ="r");
     end
  end

end

function setBC(mesh::Mesh)
  ID = zeros(Int,dof,size(mesh.coords)[2])
  IN = [] #unknown length
  cnt = 1
  #Set Dirichlet and Neumann nodes
  for i = 1:size(mesh.coords)[2]
    for j = 1:dof
      if mesh.coords[2,i] == mesh.height
        IN = [IN; cnt]
      end
      if mesh.coords[1,i] < 10e-08 #set left edge dofs Dirichlet
        ID[j,i] = 0
      else
        ID[j,i] = cnt
        cnt+=1
      end
    end #end for j
  end #end for i

  nelem = size(mesh.elems)[2]
  nnode = 4
  LM = zeros(Int, dof*nnode, nelem)
  for i = 1:nelem
    for j = 1:dof
      for k = 1:nnode
        LM[(k-1)*dof+j,i] = ID[j,mesh.elems[k,i]]
      end
    end
  end

  return MeshWithBC(mesh, ID, LM, IN);
end


function plotDisplacements(mBC::MeshWithBC, x)
  m = mBC.mesh
  #Plot original mesh
  for i = 1:size(m.coords)[2]
    plot(m.coords[1,i],m.coords[2,i],marker="x",color ="grey")
  end
  for i = 1:size(m.elems)[2]
    plot(m.coords[1,m.elems[1:4,i]]',m.coords[2,m.elems[1:4,i]]',color = "grey", linestyle="--")
  end
  #Plot displaced mesh
  cmove = zeros(size(m.coords))
  for i = 1:size(cmove)[2]
    cmove[1,i] = m.coords[1,i]
    if mBC.ID[1,i] > 0
      cmove[1,i] += x[mBC.ID[1,i]]
    end
    cmove[2,i] = m.coords[2,i]
    if mBC.ID[2,i] > 0
      cmove[2,i] += x[mBC.ID[2,i]]
    end
  end
  for i = 1:size(m.coords)[2]
    plot(cmove[1,i],cmove[2,i],marker="o",color ="b")
  end
  for i = 1:size(m.elems)[2]
    plot(cmove[1,m.elems[1:4,i]]',cmove[2,m.elems[1:4,i]]',color = "b", linestyle="-")
  end
end



end  #End Mesh

