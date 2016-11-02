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

    def euler(self):
        f = self.equation.f
        x = [self.x0]
        ym = [self.y0]
        xn = x[0]
        yn = ym[0]

        r = ode(f=f,jac=None).set_integrator('dop853',method='bdf')
        r.set_initial_value(self.y0,self.x0)
        r.integrate(xn)
        y = [r.y[0]]
        print('%(xn).5f\t%(yn).5f\t%(ya).5f\t%(e).5f'%{'xn':xn, 'yn':yn, 'ya':r.y,'e':(r.y-yn)*100/r.y})

        while 1:
            yn = yn + self.h*f(xn,yn)
            xn = xn + self.h
            r.integrate(xn)
            print('%(xn).5f\t%(yn).5f\t%(ya).5f\t%(e).5f'%{'xn':xn, 'yn':yn, 'ya':r.y,'e':(r.y-yn)*100/r.y})
            x.append(xn)
            ym.append(yn)
            y.append(r.y[0])
            if xn > self.xf:
                break
        self.plot(x,ym,y,'Euler','scatter')

    def midpoint(self):
        f = self.equation.f
        h = self.h
        x = [self.x0]
        ym = [self.y0]
        xn = x[0]
        yn = ym[0]
        
        r = ode(f=f,jac=None).set_integrator('dop853',method='bdf')
        r.set_initial_value(self.y0,self.x0)
        r.integrate(xn)
        y = [r.y[0]]
        print('%(xn).5f\t%(yn).5f\t%(ya).5f\t%(e).5f'%{'xn':xn, 'yn':yn, 'ya':r.y,'e':(r.y-yn)*100/r.y})
        
        while 1:
            yn = yn + h*f(xn+(h/2),yn + (h*f(xn,yn)/2))
            xn = xn + h
            
            r.integrate(xn)
            print('%(xn).5f\t%(yn).5f\t%(ya).5f\t%(e).5f'%{'xn':xn, 'yn':yn, 'ya':r.y,'e':(r.y-yn)*100/r.y})
        
            x.append(xn)
            ym.append(yn)
            y.append(r.y[0])
            if xn > self.xf:
                break
        self.plot(x,ym,y,'Midpoint','scatter')

    def rk4(self):
        f = self.equation.f
        h = self.h
        x = [self.x0]
        ym = [self.y0]
        xn = x[0]
        yn = ym[0]
        xf = self.xf

        r = ode(f=f,jac=None).set_integrator('dop853',method='bdf')
        r.set_initial_value(self.y0,self.x0)
        r.integrate(xn)
        y = [r.y[0]]
        print('%(xn).5f\t%(yn).5f\t%(ya).5f\t%(e).5f'%{'xn':xn, 'yn':yn, 'ya':r.y,'e':(r.y-yn)*100/r.y})
                
        while xn < xf:
            k1 = h*f(xn,yn)
            k2 = h*f(xn + h/2,yn + k1/2)
            k3 = h*f(xn + h/2,yn + k2/2)
            k4 = h*f(xn + h,yn + k3)
            xn = xn + h
            yn = yn + k1/6 + k2/3 + k3/3 + k4/6
            r.integrate(xn)
            print('%(xn).5f\t%(yn).5f\t%(ya).5f\t%(e).5f'%{'xn':xn, 'yn':yn, 'ya':r.y,'e':(r.y-yn)*100/r.y})
            x.append(xn)
            ym.append(yn)
            y.append(r.y[0])
        self.plot(x,ym, y,'4th Order RK','scatter')

    def rk45(self):
        f = self.equation.f
        h = self.h
        x = [self.x0]
        ym = [self.y0]
         
        xn = x[0]
        yn = ym[0]
        xf = self.xf
        tol = self.tol
        alpha = self.alpha
        
        r = ode(f=f,jac=None).set_integrator('dop853',method='bdf')
        r.set_initial_value(self.y0,self.x0)
        r.integrate(xn)
        y = [r.y[0]]
        print('%(xn).5f\t%(yn).5f\t%(ya).5f\t%(e).5f'%{'xn':xn, 'yn':yn, 'ya':r.y,'e':(r.y-yn)*100/r.y})
        
        while xn < xf:
            h = min(h, xf - xn, self.hmax)

            k1 = h*f(xn,yn)
            k2 = h*f(xn + h/4, yn + k1/4)
            k3 = h*f(xn + 3*h/8, yn + 3*k1/32 + 9*k2/32)
            k4 = h*f(xn + 12*h/13, yn + 1932*k1/2197 - 7200*k2/2197 + 7296*k3/2197 )
            k5 = h*f(xn + h, yn + 439*k1/216 - 8*k2 + 3680*k3/513 - 845*k4/4104)
            k6 = h*f(xn + h/2, yn - 8*k1/27 + 2*k2 - 3544*k3/2565 + 1859*k4/4104 - 11*k5/40)

            y4 = yn + 25*k1/216 + 1408*k3/2565 + 2197*k4/4101 - k5/5
            y5 = yn + 16*k1/135 + 6656*k3/12825 + 28561*k4/56430 - 9*k5/50 + 2*k6/55

            # k1 = h*f(xn,yn)
            # k2 = h*f(xn + h/2,yn + k1/2)
            # k3 = h*f(xn + h/2,yn + k2/2)
            # k4 = h*f(xn + h,yn + k3)
            # y4 = yn + k1/6 + k2/3 + k3/3 + k4/6

            # k1 = h*f(xn, yn)
            # k2 = h*f(xn + h/4, yn + k1/4)
            # k3 = h*f(xn + h/4, yn + k1/8 + k2/8)
            # k4 = h*f(xn + h/2, yn - k2/2 + k3 )
            # k5 = h*f(xn + 3*h/4, yn + 3*k1/16 + 9*k4/16)
            # k6 = h*f(xn + h, yn - 3*k1/7 + 2*k2/7 + 12*k3/7 - 12*k4/7 + 8*k5/7)
            # y5 = yn + 7*k1/90 + 32*k3/90 + 12*k4/90 + 32*k5/90 + 7*k6/90

            

            E = abs(y5 - y4)
            # print(R, abs(zn-yn), h)
            
            if E < tol:
                yn = y5
                xn = xn + h
                r.integrate(xn)
                print('%(xn).5f\t%(yn).5f\t%(ya).5f\t%(e).5f'%{'xn':xn, 'yn':yn, 'ya':r.y,'e':(r.y-yn)*100/r.y})
                x.append(xn)
                ym.append(yn)
                y.append(r.y[0])
            else:
                h = 0.84*h*(tol*h/E)**alpha 
                
        self.plot(x,ym,y,"RK45","scatter")

    def plot(self,x,ym,y,method,pl):
        import matplotlib.pyplot as plt
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_title(method+" Method")
        ax.set_ylabel("Y")
        ax.set_xlabel("X")
        ax.plot(x,y,'-',label='Analytic')
        ax.scatter(x,ym, s=30, c='r', marker="s", label=method)
        plt.legend(loc='best')
        plt.show()
        fig.savefig(method+" plot.png")
    
o = ODESolver()
o.euler()
o.midpoint()
o.rk4()
o.rk45()