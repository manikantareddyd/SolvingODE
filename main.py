import sys
import parser
from math import *
from sympy import *
import numpy as np


from equation import *

class ODESolver:
    def __init__(self):
        self.equation = Equation()
        self.x0 = float(input())
        self.y0 = float(input())
        self.xf = float(input())
        self.h = float(input())
        self.hmax = float(input())
        self.alpha = float(input())
        self.tol = float(input())

    def euler(self):
        x = [self.x0]
        y = [self.y0]
        xn = x[0]
        yn = y[0]
        print('%(xn).4f\t%(yn).4f'%{'xn':xn, 'yn':yn})
        while 1:
            yn = yn + self.h*self.equation.f(xn,yn)
            xn = xn + self.h
            print('%(xn).4f\t%(yn).4f'%{'xn':xn, 'yn':yn})
            x.append(xn)
            y.append(yn)
            if xn > self.xf:
                break
        
o = ODESolver()
o.euler()