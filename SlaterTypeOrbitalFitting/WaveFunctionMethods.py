from numpy import pi, exp
from math import factorial


class SphericalHarmonic:
    def __init__(self,l,m):
        self.l = l
        self.m = m
    def Y_lm(self, x, y, z):
        r = ((x**2)+(y**2)+(z**2))**0.5
        if r==0:
            return 0
        else:
            if self.l==0 and self.m==0:
                return (1/2)*((1/pi**0.5))
            elif self.l==1 and self.m==-1:
                return ((3/(4*pi))**0.5)*(y)/r
            elif self.l==1 and self.m==0:
                return ((3/(4*pi))**0.5)*(z)/r
            elif self.l==1 and self.m==1:
                return ((3/(4*pi))**0.5)*(x)/r
            elif self.l==2 and self.m==-2:
                return (1/2)*((15/pi)**0.5)*(y)*(x)/(r**2)
            elif self.l==2 and self.m==-1:
                return (1/2)*((15/pi)**0.5)*(y)*(z)/(r**2)
            elif self.l==2 and self.m==0:
                return (1/4)*((5/pi)**0.5)*(3*(z)**2-r**2)/(r**2)
            elif self.l==2 and self.m==1:
                return (1/2)*((15/pi)**0.5)*(z)*(x)/(r**2)
            elif self.l==2 and self.m==2:
                return (1/4)*((15/pi)**0.5)*((x)**2-(y)**2)/(r**2)
            else:
                return 0


class RadialPart:
    def __init__(self,n,zeta):
        self.N=((2*zeta)**n)*((2*zeta/(factorial((2*n))))**0.5)
        self.n=n
        self.Zeta=zeta
    def CalcRad(self,x,y,z):
        R=(((x)**2)+((y)**2)+((z)**2))**0.5
        return self.N*(R**(self.n-1))*(exp(-self.Zeta*R))
    
    
class SlaterTypeOrbital:
    def __init__(self,l,m,n,zeta):
        self.Radial = RadialPart(n,zeta)
        self.Harmonic = SphericalHarmonic(l,m)
    def CalcWaveFunc(self, x, y, z):
        RValue = self.Radial.CalcRad(x, y, z)
        YValue = self.Harmonic.Y_lm(x, y, z)
        return RValue*YValue

