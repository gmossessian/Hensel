'''
Created on Jul 9, 2015

@author: George Mossessian

Hensel's Lemma: 
    let f(x) be a polynomial with integer coefficients, p a prime number, and k>=2, 
    and consider equation f(x) mod p^k = 0. Suppose r is a solution to f(x) mod p^{k-1} = 0. Then:
    I) if f'(r) mod p != 0, then there is a unique integer t, with 0<=t<p, such that 
        f(r+tp^(k-1)) mod p^k = 0, given by t=-(f'(r)^(-1))(f(r)/p^(k-1)) mod p.
    II) if f'(r) mod p = 0 and f(r) mod p^k = 0, then f(r+tp^(k-1)) mod p^k = 0 for all integers t.
    III) if f'(r) mod p = 0 and f(r) mod p^k != 0, then f(x) mod p^k = 0 has no solutions with x = r mod p^(k-1).
    
INPUT: int p, int k, a_n, a_(n-1), ... , a_0 
    where the congruence equation in question is 
    a_nx^n+a_(n-1)x^(n-1)+...+a_0x^0 (mod p^k) = 0.
    
    so to solve x^2 + x + 7 = 0 mod 25, you would type on the command-line
    
    Hensel 5 2 1 1 7
'''

from sys import argv

class Polynomial(list):
    """A polynomial whose coeffients, from highest power to lowest, are given by list."""
    def __init__(self, clist):
        self.coeffs = clist
    
    def __repr__(self):
        sign = {
            (True, True): '-',
            (True, False): '',
            (False, True): ' - ',
            (False, False): ' + '
        }
        
        poly = []
        
        for n, a in reversed(list(enumerate(reversed(self.coeffs)))):
            s=sign[not poly, a<0]
            a=abs(a)
            if a==1 and n !=0:
                a = ''
            f={0:'{}{}', 1: '{}{}x'}.get(n, '{}{}x^{}')
            if a != 0: poly.append(f.format(s,a,n))
        return ''.join(poly) or ''
        
def FormalDerivative(f):
    """Returns a Polynomial which is the formal derivative of the polynomial f"""
    dc = Polynomial([])
    l=len(coeffs)
    for i in range(l-1):
        dc.coeffs.append((l-i-1)*coeffs[i])
    return dc

def Evaluate(f, x):
    """Evaluate f(x)"""
    r=0
    n=len(f.coeffs)-1
    for i in range(n+1):
        r += f.coeffs[i]*(x**(n-i))
    return r

def ExtendedEuclideanAlgorithm(a,b,c):
    """Returns the least x, y such that ax + by = c."""
    revflag=False
    if b>a:
        a,b=b,a
        revflag=True
    r , q = [a,b,a%b] , int(a/b)
    s , t = [1,0,1] , [0,1,-q]
    if r[2]==0:
        if c%b == 0:
            if revflag == False: return 0,c/b
            return c/b, 0
        return False
    while r[2]>0:
        q=int(r[1]/r[2])
        r[0], r[1], r[2] = r[1], r[2], r[1]%r[2]
        s[0], s[1], s[2] = s[1], s[2], s[1]-s[2]*q
        t[0], t[1], t[2] = t[1], t[2], t[1]-t[2]*q
    if c%r[1] != 0:
        return False
    d=c/r[1]
    if revflag: 
        s,t = t,s
    return d*s[1],d*t[1]

def Hensel(f,p,k):
    """Returns solutions to f(x) mod p^k = 0. """
    if k < 1: return False
    if k == 1:
        r=[]
        for i in range(p):
            if Evaluate(f,i)%p == 0:
                r.append(i)
        return r
    r = Hensel(f,p,k-1)
    df=FormalDerivative(f)
    rnew=[]
    for i,n in enumerate(r):
        dfr=Evaluate(df,n)
        fr=Evaluate(f,n)
        if dfr%p != 0:
            t=(-(ExtendedEuclideanAlgorithm(dfr, p, 1)[0])*int(fr/p**(k-1)))%p
            rnew.append(r[i]+t*p**(k-1))
        if dfr%p == 0:
            if fr % p**k == 0:
                for t in range(0,p):
                    rnew.append(r[i]+t*p**(k-1))                
    return rnew

if __name__ == '__main__':
    p = int(argv[1])
    k = int(argv[2])
    coeffs = []
    for i in range(3,len(argv)):
        coeffs.append(int(argv[i]))
    f = Polynomial(coeffs)
    print('Using Hensel\'s Lemma to find solutions for: ')
    print(f, ' mod ', p, '^' ,k,'= 0', sep="")
    r = Hensel(f,p,k)
    print('Solutions:',r)