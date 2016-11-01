from scipy.integrate import ode
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
        self.x = [self.x0]
        self.y = [self.y0]

    def analytic(self):
        foo = self.equation.f
        r = ode(f=foo,jac=None).set_integrator('dop853',method='bdf')
        r.set_initial_value(self.y0,self.x0)
        while r.successful() and r.t < self.xf:
            r.integrate(r.t + 0.05)
            self.x.append(r.t)
            self.y.append(r.y)
        self.plot(self.x,self.y,'Analytic ','plot')

    def euler(self):
        x = [self.x0]
        y = [self.y0]
        xn = x[0]
        yn = y[0]
        print('%(xn).5f\t%(yn).5f'%{'xn':xn, 'yn':yn})
        while 1:
            yn = yn + self.h*self.equation.f(xn,yn)
            xn = xn + self.h
            print('%(xn).5f\t%(yn).5f'%{'xn':xn, 'yn':yn})
            x.append(xn)
            y.append(yn)
            if xn > self.xf:
                break
        self.plot(x,y,'Euler','scatter')

    def midpoint(self):
        f = self.equation.f
        h = self.h
        x = [self.x0]
        y = [self.y0]
        xn = x[0]
        yn = y[0]
        print('%(xn).5f\t%(yn).5f'%{'xn':xn, 'yn':yn})
        while 1:
            yn = yn + h*f(xn+(h/2),yn + (h*f(xn,yn)/2))
            xn = xn + h
            print('%(xn).5f\t%(yn).5f'%{'xn':xn, 'yn':yn})
            x.append(xn)
            y.append(yn)
            if xn > self.xf:
                break
        self.plot(x,y,'Midpoint','scatter')

    def rk4(self):
        f = self.equation.f
        h = self.h
        x = [self.x0]
        y = [self.y0]
        xn = x[0]
        yn = y[0]
        print('%(xn).5f\t%(yn).5f'%{'xn':xn, 'yn':yn})
        while 1:
            xm = xn + h/2
            k1 = f(xn,yn)
            ym = yn + h*k1/2
            k2 = f(xm,ym)
            ym1 = yn + h*k2/2
            k3 = f(xm,ym1)
            ye = yn + h*k3
            xe = xn + h
            k4 = f(xe,ye)
            phi = (k1 + 2*k2 + 2*k3 + k4)/6

            xn = xe
            yn = yn + h*phi
            print('%(xn).5f\t%(yn).5f'%{'xn':xn, 'yn':yn})
            x.append(xn)
            y.append(yn)
            if xn > self.xf:
                break
        self.plot(x,y,'4th Order RK','scatter')

    def plot(self,x,y,method,pl):
        import matplotlib.pyplot as plt
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_title(method+" Method")
        ax.set_ylabel("Y")
        ax.set_xlabel("X")
        if pl == 'plot':
            ax.plot(x,y,'-')
        else:
            ax.scatter(x,y, s=30, c='r', marker="s", label=method)
        plt.legend(loc='best')
        plt.show()
        fig.savefig(method+" plot.png")
    
o = ODESolver()
o.euler()
o.midpoint()
o.rk4()