module quadrature

export integrate,gauss

type QuadRule
  x
  w
  
  function QuadRule(x,w)
    new(x,w)
  end
end

#Construct one-dimensional rules
function gauss(n::Int)
    x = zeros(n)
    w = zeros(n)

    m = floor((n+1)/2) 
    xm = 0.0;    xl = 1.0
    
    for i = 1:m
        z  = cos(pi*(i-0.25)/(n+0.5))
        pp = 0
        while true
            p1 = 1.0 ;    p2 = 0.0
            for j = 1:n
                p3 = p2;  p2 = p1
                p1 = ((2.0*j-1.0)*z*p2-(j-1.)*p3)/j
            end
            pp = n*(z*p1-p2)/(z*z-1.0)
            z1 = z
            z  = z1-p1/pp
            if abs(z-z1)<eps()
              break
            end
        end
        x[i] = xm-xl*z
        x[n+1-i] = xm+xl*z
        w[i] = 2.0*xl/((1.0-z*z)*pp*pp)
        w[n+1-i] = w[i]
    end
    #Transform to [0,1]
    w = w/2.;
    x = (x+1.)/2.

    q = QuadRule(x,w)
    return q
end

#Construct 2d rule on the cube [0,1]^2 given two one dimensional rules
function tensor(q1::QuadRule, q2::QuadRule)
  l1 = size(q1.x)[1]
  l2 = size(q2.x)[1]
  z = zeros(l1*l2,2)
  w = zeros(l1*l2)

  #slow changing q1.x
  z[:,1] = repmat(q1.x,l2,1)
  w = repmat(q2.w, l2, 1)
  #fast changing q2.x
  y1 = repmat(q2.x',l1,1)
  z[:,2] = reshape(y1,l1*l2,1)
  v1 = repmat(q2.w',l1,1)
  w = reshape(v1, l1*l2, 1).*w

  return QuadRule(z,w)
end

#Transform the quadrature rule to the element
function transform2Element(q::QuadRule,element)
  A = hcat(element[:,2]-element[:,1], element[:,3]-element[:,1])
  y = zeros(size(q.x))
  for i = 1:size(q.x)[1]
    y[i,:] = (element[:,1] + A*q.x[i,:]')'
  end
  return QuadRule(y,q.w*abs(det(A)))
end

#Integrate the function f on the element
function integrate(f, element,intorder)
  q = gauss(intorder)
  q2 = tensor(q,q)
  q3 = transform2Element(q2, element)
  return f(q3.x)'*q3.w
end

#Integrate the function f on the reference element
function integrate(f,intorder)
  q = gauss(intorder)
  q2 = tensor(q,q)
  help_f = f(q2.x)
  if size(help_f) == size(q2.w)
    return help_f'*q2.w
  elseif size(help_f)[1] == size(q2.w)[1]
    fh = reshape(help_f,size(q2.w))
    return fh'*q2.w
  else
    print("Wrong input size to quadrature \n")
  end
end

end
