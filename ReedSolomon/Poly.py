from GF import GF256
class Polynomial(object):#highest at front
    def __init__(self, coef = ()):
        c=list(coef)
        while c and c[0]==0:
            c.pop(0)
        if not c:
            c.append(0)
        self.coef = tuple(c)
    def __len__(self):
        return len(self.coef)
    def degree(self):
        return len(self.coef)-1
    def __add__(self, other):
        diff = len(self) - len(other)
        if diff > 0:
            t1 = self.coef
            t2 = (0,) * diff + other.coef
        else:
            t1 = (0,) * (-diff) + self.coef
            t2 = other.coef
        return self.__class__(x+y for x,y in zip(t1, t2))
    def __neg__(self):
        return self.__class__(-a for a in self.coef)
    def __sub__(self,other):
        return self+(-other)
    def __mul__(self, other):
        terms = [0] * (len(self) + len(other))
        for i1, c1 in enumerate(reversed(self.coef)):
            if c1 == 0:
                continue
            for i2, c2 in enumerate(reversed(other.coef)):
                terms[i1+i2] += c1*c2

        return self.__class__(reversed(terms))
    
    def __divmod__(dividend, divisor):#low efficiency
        class_ = dividend.__class__
        dividend_power = dividend.degree()
        dividend_coefficient = dividend.coef[0]

        divisor_power = divisor.degree()
        divisor_coefficient = divisor.coef[0]

        quotient_power = dividend_power - divisor_power
        if quotient_power < 0:

            return class_((0,)), dividend

        quotient_coefficient = dividend_coefficient // divisor_coefficient
        quotient = class_( (quotient_coefficient,) + (0,) * quotient_power )

        remander = dividend - quotient * divisor

        if remander.coef == (0,):
            return quotient, remander


        morequotient, remander = divmod(remander, divisor)
        return quotient + morequotient, remander
    def __floordiv__(self,other):
        return divmod(self,other)[0]
    def __mod__(self,other):
        return divmod(self,other)[1]
    def __eq__(self,other):
        return self.coef == other.coef
    def __ne__(self,other):
        return self.coef != other.coef
    

    def get_coefficient(self, degree):
        if degree > self.degree():
            return 0
        else:
            return self.coef[-(degree+1)]
    def coefficient(self):
        return self.coef
    def evaluate(self, x):
        c =0
        p =1
        for term in reversed(self.coef):
            c = c + term * p
            p = p * x

        return c


