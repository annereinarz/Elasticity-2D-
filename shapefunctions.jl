module shapefunctions

export N, dN

# a node
function N(a,x,y)
  if  a == 1
    return (1-x).*(1-y)
  end
  if a == 2
    return x.*(1-y)
  end
  if a == 3
    return x.*y
  end
  if a == 4
    return (1-x).*y
  end
  throw("Squares have only four corners")
end

#i derivative direction
function dN(a::Int,x,y,i::Int)
  if a == 1
    if i == 1
      return -(1-y)
    end
    if  i == 2
      return -(1-x)
    end
  end
  if a == 2
    if i == 1
      return (1-y)
    end
    if  i == 2
      return -x
    end
  end
  if a == 3
    if i == 1
      return y
    end
    if  i == 2
      return x
    end
  end
  if a == 4
    if i == 1
      return -y
    end
    if  i == 2
      return (1-x)
    end
  end
  throw(string("Squares have only four corners", a))
end #end dN 

end #end module
