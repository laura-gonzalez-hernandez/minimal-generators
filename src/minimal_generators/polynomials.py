class Monomial:
    def __init__(self, coefficient:int, exponents:list[int], variables=['x','y','z']):
        self.sign='+' if coefficient>=0 else '-'
        self.coefficient=abs(coefficient)
        self.variables=variables
        self.exponents=exponents
    def write(self)->str:
        out=''
        if self.coefficient not in [0,1]:
            out+=str(self.coefficient)
        for v,e in zip(self.variables,self.exponents):
            if e!=0:
                out+=v+str(e if e>1 else '')
        return out
    def write_singular(self)->str:
        out=''
        if self.coefficient not in [0,1]:
            out+=str(self.coefficient)
        for v,e in zip(self.variables,self.exponents):
            if e!=0:
                out = out + ('' if out=='' else '*') + v+'^'+str(e)
        return out

class Polynomial:
    def __init__(self, monomials:list[Monomial]):
        self.monomials=monomials
    def write(self, singular=False)->str:
        out=''
        for m in self.monomials:
            if m.coefficient!=0:
                mono=m.write_singular() if singular else m.write()
                out += (mono+' ' if out=='' and m.sign=='+' else m.sign+mono+' ')
        return out + (';' if singular else '')
    def print(self, singular=False, name=None):
        txt=self.write(singular=singular)
        print(f'{name} = {txt}' if name else txt)
    def __add__(self, other:'Polynomial')->'Polynomial':
        return Polynomial(self.monomials+other.monomials)

def build_monomials(coef:list[int], exp_x:list[int], exp_y:list[int], exp_z:list[int]) -> list[Polynomial]:  # construeix objectes Monomial/Polynomial a partir de coef i exponents
    """
    Build list of Polynomial objects from coefficient and exponent lists.
    """
    polys=[]
    for c, ex, ey, ez in zip(coef, exp_x, exp_y, exp_z):
        mono=[]
        for coeff,x,y,z in zip(c, ex, ey, ez):
            mono.append(Monomial(coeff,[x,y,z]) if coeff!=0 else Monomial(0,[0,0,0]))
        polys.append(Polynomial(mono))
    return polys

def build_f_polys(n:int, g:list[Polynomial], h1:list[Polynomial], h2:list[Polynomial], h3:list[Polynomial], h4:list[Polynomial]) -> list[Polynomial]:
    """
    Build the list of f polynomials from g and h polynomials.
    If n is even, insert specific polynomials at index 2."""
    f=[]
    if n%2==0:
        g.insert(n-3, Polynomial([Monomial(1,[1,n-3,1]), Monomial(-1,[0,n-1,0]), Monomial(-1,[0,0,n-1])]))
        h1.insert(n-3, Polynomial([Monomial(0,[0,0,0])]))
        h2.insert(n-3, Polynomial([Monomial(0,[0,0,0])]))
        h3.insert(n-3, Polynomial([Monomial(0,[0,0,0])]))
        h4.insert(n-3, Polynomial([Monomial(0,[0,0,0])]))
    for idx in range(len(g)):
        f.append(g[idx]+h1[idx]+h2[idx]+h3[idx]+h4[idx])
    return f