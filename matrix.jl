dof = 2

function globalMatrix(mD,intorder)
  #loop over elements, create local stiffness matrix
  #and place into global matrix
  K   = zeros(maximum(mD.ID),maximum(mD.ID))
  for i = 1:size(mD.mesh.elems)[2]
    element = mD.mesh.coords[:,mD.mesh.elems[:,i]]
    elemdof = size(mD.LM)[1]
    KPQ = localHourglass(element,mD.mesh.height,intorder)
    #equation numbers, write appropriate matrix entries
    ih = [1:elemdof]
    local_ind = ih[mD.LM[:,i] .!= 0]
    global_ind = mD.LM[local_ind,i]
    K[global_ind,global_ind] += KPQ[local_ind,local_ind]
  end
  return K
end


using quadrature
function localHourglass(element,height,intorder)
  a = 1./5
  Kelem = zeros(8,8)
  Kelem1 = zeros(8,8)
  Kelem2 = zeros(8,8)
  A = hcat(element[:,2]-element[:,1], element[:,3]-element[:,1])
  elem_center = element[:,1] + A*[0.5, 0.5]
  elem_det = abs(det(A))

  for a = 1:4   #run over all local nodes
    for i = 1:dof
      for b = 1:4
        for j = 1:dof
          integral = quadrature.integrate(x->f(a,b,i,j,x,elem_center,height),intorder)
          Kelem1[dof*(a-1)+i, dof*(b-1)+j] = integral[1]*elem_det

          integral = quadrature.integrate(x->f(a,b,i,j,x,elem_center,height),intorder+1)
          Kelem2[dof*(a-1)+i, dof*(b-1)+j] = integral[1]*elem_det
        end
      end
    end
  end
  
  lambda2,v2 = eig(Kelem2); lambda2 = real(lambda2); v2 = real(v2)
  lambda1 = diag(v2'*Kelem1*v2)
  #l = setdiff(round(lambda2,8),round(lambda1,8))
  for i = 1:8
    if i==6 || i==5  #spurious modes
      lambda1[i] = lambda2[i]*a
    end
    Kelem += v2[:,i]*lambda1[i]*v2[:,i]'
  end
  return Kelem
end

function localMatrix(element,height,intorder)
  Kelem = zeros(8,8)
  A = hcat(element[:,2]-element[:,1], element[:,3]-element[:,1])
  elem_center = element[:,1] + A*[0.5, 0.5]
  elem_det = abs(det(A))
  for a = 1:4   #run over all local nodes
    for i = 1:dof
      for b = 1:4
        for j = 1:dof
           integral = quadrature.integrate(x->f(a,b,i,j,x,elem_center,height),intorder)
           Kelem[dof*(a-1)+i, dof*(b-1)+j] = integral[1]*elem_det
        end
      end
    end
  end
  return Kelem
end

using shapefunctions
function f(a::Int,b::Int,i,j,x,elem_center,height)
  sol = zeros(size(x)[1],1); 
  D = Cijkl(x,elem_center,height)
  for c = 1:size(x)[1]
    Ba = [shapefunctions.dN(a,x[c,1],x[c,2],1) 0;
          0                                    shapefunctions.dN(a,x[c,1],x[c,2],2);
          shapefunctions.dN(a,x[c,1],x[c,2],2) shapefunctions.dN(a,x[c,1],x[c,2],1)]
    Bb = [shapefunctions.dN(b,x[c,1],x[c,2],1) 0;
          0                                    shapefunctions.dN(b,x[c,1],x[c,2],2);
          shapefunctions.dN(b,x[c,1],x[c,2],2) shapefunctions.dN(b,x[c,1],x[c,2],1)]
    sol[c] = (Ba'*D*Bb)[i,j]
  end
  return sol;
end

function Cijkl(x,elem_center,height)
  #set up material properties
  num_layers = 5
  mu = 4.083;  lambda  = 1.38
  D = [lambda+2*mu lambda 0; lambda (lambda+2*mu)/10 0; 0 0 mu ]
  for i = 2:2:num_layers
    if (i-1)/num_layers*height < elem_center[2] && elem_center[2] < i/num_layers*height      
      #Youngs modulus E = .5, Poisson ratio nu = 0.2
      #lambda = E*nu/((1+nu)*(1-2*nu)) = .1/(1.2*0.6) = 1.38/10
      #mu = E/(2+2*nu) = .5/(2+2*0.2) = .5/2.4 = 2.083/10 ,
      mu = 2.083/10;  lambda = 1.38/10
      D = [lambda+2*mu lambda      0;
      lambda      lambda+2*mu 0;
      0           0           mu]
    end
  end
  return D 
end
