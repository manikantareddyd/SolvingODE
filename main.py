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
        xf = self.xf
        print('%(xn).5f\t%(yn).5f'%{'xn':xn, 'yn':yn})
        while xn < xf:
            k1 = h*f(xn,yn)
            k2 = h*f(xn + h/2,yn + k1/2)
            k3 = h*f(xn + h/2,yn + k2/2)
            k4 = h*f(xn + h,yn + k3)
            xn = xn + h
            yn = yn + k1/6 + k2/3 + k3/3 + k4/6
            print('%(xn).5f\t%(yn).5f'%{'xn':xn, 'yn':yn})
            x.append(xn)
            y.append(yn)
        self.plot(x,y,'4th Order RK','scatter')

    def rk45(self):
        f = self.equation.f
        h = self.h
        x = [self.x0]
        y = [self.y0]
        xn = x[0]
        yn = y[0]
        xf = self.xf
        tol = self.tol
        alpha = self.alpha
        print('%(xn).5f\t%(yn).5f'%{'xn':xn, 'yn':yn})
        while xn < xf:
            h = min(h, xf - xn, self.hmax)
            # print('t',xn,h)
            k1 = h*f(xn,yn)
            k2 = h*f(xn + h/4, yn + k1/4)
            k3 = h*f(xn + 3*h/8, yn + 3*k1/32 + 9*k2/32)
            k4 = h*f(xn + 12*h/13, yn + 1932*k1/2197 - 7200*k2/2197 + 7296*k3/2197 )
            k5 = h*f(xn + h, yn + 439*k1/216 - 8*k2 + 3680*k3/513 - 845*k4/4104)
            k6 = h*f(xn + h/2, yn - 8*k1/27 + 2*k2 - 3544*k3/2565 + 1859*k4/4104 - 11*k5/40)

            xn = xn + h
            y4 = yn + 25*k1/216 + 1408*k3/2565 + 2197*k4/4101 - k5/5
            y5 = yn + 16*k1/135 + 6656*k3/12825 + 28561*k4/56430 - 9*k5/50 + 2*k6/55

            E = abs(y5 - y4)
            # print(R, abs(zn-yn), h)
            
            if E < tol:
                yn = y5
                print('%(xn).5f\t%(yn).5f'%{'xn':xn, 'yn':yn})
                x.append(xn)
                y.append(yn)
            else:
                delta = (tol/E)**alpha
                h = delta*h
                
        self.plot(x,y,"RK45","scatter")

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
# o.euler()
# o.midpoint()
o.rk4()
# o.rk45()