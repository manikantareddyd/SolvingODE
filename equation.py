import sys
import parser
from math import *
from sympy import *
import numpy as np

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

    def f(self,x0,y0):
        x = Symbol('x')
        y = Symbol('y')
        z = sympify(self.equation)
        z = z.subs(x,x0*1.0)*1.0
        val = z.subs(y,y0*1.0)*1.0
        return val

    # def df(self,voo):
    #     x = Symbol('x')
    #     y = sympify(self.equation)
    #     yprime = diff(y,x)
    #     # f = lambdify(x, yprime, 'numpy')
    #     val = yprime.subs(x,voo*1.0)*1.0
    #     return val