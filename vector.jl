using shapefunctions
dof = 2

function globalVector(mD,f, h,intorder)
  #loop over elements, create local stiffness matrix
  #and place into global matrix
  B   = zeros(maximum(mD.ID))
  elemdof = size(mD.LM)[1]
  for i = 1:size(mD.mesh.elems)[2]
    #Get element coordinates
    element = mD.mesh.coords[:,mD.mesh.elems[:,i]]
    BP = elementVector(f,element,intorder)
    BNeu = neumannVector(h,element,mD)
    #Set index sets
    ih = [1:elemdof]
    local_ind = ih[mD.LM[:,i] .!= 0]
    global_ind = mD.LM[local_ind,i]
    #equation numbers, write appropriate vector entries
    B[global_ind] += BP[local_ind] +BNeu[local_ind] 
  end
  return B
end

function neumannVector(h,element,mD)
  Bneu = zeros(8)
  for a = 1:4
    for i = 1:dof
      if abs(element[1,a] - mD.mesh.width) < 1e-06 && i == 1
      #if abs(element[2,a] - mD.mesh.height) < 1e-06 && i == 1
        Bneu[dof*(a-1)+i] = h
      end
    end
  end
 return Bneu
end

function elementVector(f,element,intorder)
  A = hcat(element[:,2]-element[:,1], element[:,3]-element[:,1])
  b = element[:,1]
  Belem = zeros(8)
  #right hand side
  for a = 1:4  #loop over all elements
    for i = 1:dof
      value = integrate(x->f(x).*reshape(shapefunctions.N(a, x[:,1],x[:,2]),size(f(x))),intorder)
      Belem[dof*(a-1)+i] +=  value[1,1]*abs(det(A))
    end
  end
  return Belem
end

