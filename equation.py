import sys
from math import *
import sympy

class Equation:
    def __init__(self):
        while 1:
            try:
                print("Please enter equation using standard system notation.\nNote:The Variables should be x and y only")
                self.equation = str(input())
                # self.equation = "5*(exp(-100*(x-2)^2)) - 0.5*y"
                # self.equation = "(600*(x^4))-(500*(x^3))+(200*(x^2))-(20*x)-1"
                # self.equation = "(x^3)+(x^2)-(4*x)-4"
                if(self.equation == "exit"):
                    break
                # print(self.equation)
                break
            except:
                print("\nOops!",sys.exc_info()[0],"occured. Try again!")

    def f(self,t,y):
        x = t
        xo = sympy.Symbol('x')
        yo = sympy.Symbol('y')
        z = sympy.sympify(self.equation)
        z = z.subs(xo,x*1.0)*1.0
        val = z.subs(yo,y*1.0)*1.0
        return val
