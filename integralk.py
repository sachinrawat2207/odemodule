import numpy as np
### These methods are explicit method
class interpolation:
    def lagrangeint(x,y,k):#Functions that return the coefficient of the polynomial
        coeff=0
        for i in range(k):
            p=np.poly(x[np.arange(k)!=i])
            coeff=coeff+y[i]*p/(np.polyval(p,x[i]))
        return coeff


class derivative:
    pass
class integeration:
    pass

class ode:
    pass

def eularforward(f,x,inity,dx):
    return inity+f(x,inity)*dx

def RK2(f,x,inity,dx):
    mid=inity+f(x,inity)*dx/2
    return inity+dx*f(x+dx/2,mid)

def RK4(f,x,inity,dx):
    k1=dx*f(x,inity)
    y=inity+k1/2
    k2=dx*f(x,y)
    y=inity+k2/2
    k3=dx*f(x,y)
    y=inity+k3
    k4=dx*f(x,y)
    return inity+(k1+2*k2+2*k3+k4)/6

def pefrl(f,x,inity,dx):
    PEFRL_a=0.1786178958448091
    PEFRL_b=-0.2123418310626054
    PEFRL_c=-0.06626458266981849
    const_1=PEFRL_a*dx
    const_2=(1-2*PEFRL_b)*dx/2
    const_3=PEFRL_c*dx
    const_4=PEFRL_b*dx
    const_5=(1-2*(PEFRL_c+PEFRL_a))*dx

    #PEFRL implementation
    y=inity+const_1*f(x,inity)
    ymid=inity+const_2*f(x,y)
    y=y+const_3*f(x,ymid)
    ymid=ymid+const_4*f(x,y)
    y=y+const_5*f(x,ymid)
    ymid=ymid+const_4*f(x,y)
    y=y+const_3*f(x,ymid)
    ymid=ymid+const_2*f(x,y)
    return y+(const_1)*f(x,ymid)