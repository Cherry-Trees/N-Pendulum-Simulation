import numpy as np
from scipy.integrate import odeint


class PendulumSystem:

    g = 9.81
    class _Pendulum:

        def __init__(self, m, L, th0, thdot0=0, fric=0, experienceGravity=True,
                     showbob=True, bobcolor="DEFAULT", rodcolor="black", rodsize=2.5,
                     visualstyle="trail", visualstylecolors="DEFAULT") -> None:
            self.m = m
            self.L = L
            self.th = th0
            self.thdot = thdot0
            self.fric = fric
            self.experienceGravity = experienceGravity
            self.visualstyle = visualstyle
            self.visualstylecolors = visualstylecolors
            self.showbob = showbob

            self.trail = None
            self.tracer = None
            self.v_vector = None
            self.a_vector = None
            self.rod_sprite = None
            self.bob_sprite = None
            self.bobcolor = bobcolor
            self.rodcolor = rodcolor
            self.rodsize = rodsize

            self.th_data = []
            self.thdot_data = []
            self.thddot_data = []
            self.x_data = []
            self.y_data= []
            self.xdot_data = []
            self.ydot_data = []
            self.xddot_data = []
            self.yddot_data = []
            self.T_data = []
            self.V_data = []


    def __init__(self) -> None:
        self.pl = []
        self.n = 0
        self.L_tot = 0


    def addPendulum(self, m, L, th0, **kwargs) -> None:
        '''Add a Pendulum object to the system'''
        defaultKwargs = {
            "thdot0": 0, 
            "fric": 0,
            "experienceGravity": True,
            "visualstyle": "trail",
            "visualstylecolors": "DEFAULT",
            "showbob": True,
            "bobcolor": "DEFAULT",
            "rodcolor": "black",
            "rodsize": 2.5
        }
        kwargs = {**defaultKwargs, **kwargs}
        self.pl.append(self._Pendulum(m=m, L=L, th0=th0, thdot0=kwargs["thdot0"], fric=kwargs["fric"], 
                                      experienceGravity=kwargs["experienceGravity"], showbob=kwargs["showbob"],
                                      visualstyle=kwargs["visualstyle"], visualstylecolors=kwargs["visualstylecolors"],
                                      bobcolor=kwargs["bobcolor"], rodcolor=kwargs["rodcolor"],
                                      rodsize=kwargs["rodsize"]))
        self.n += 1
        self.L_tot += L
    

    def _ode_system(self, S, t) -> np.ndarray:

        def A():
            A = []
            for i in range(self.n):
                row = []
                for k in range(self.n):
                    mqi = 0
                    for q in range(self.n):
                        if q>=k and i<=q:
                            mqi += self.pl[q].m
                    row.append(mqi*self.pl[i].L*self.pl[k].L*np.cos(self.pl[i].th - self.pl[k].th))
                A.append(row)
            return A

        def B():
            B = []
            for i in range(self.n):
                B_i = 0
                for k in range(self.n):
                    mqi = 0
                    B_i -= (
                        self.g*self.pl[i].L*np.sin(self.pl[i].th)*self.pl[k].m*
                        int(i<=k)*int(self.pl[k].experienceGravity)
                    ) + self.pl[i].fric*self.pl[i].thdot
                    for q in range(self.n):
                        if q>=k and i<=q:
                            mqi += self.pl[q].m
                    B_i -= (
                        mqi*self.pl[i].L*self.pl[k].L*np.sin(self.pl[i].th - self.pl[k].th)*
                        (self.pl[k].thdot**2)
                    )
                B.append(B_i)
            return B

        for i in range(self.n):
            self.pl[i].th = S[i]
            self.pl[i].thdot = S[i+self.n]

        thddot_data = np.linalg.solve(A(), B()).tolist()
        thdot_data = [p.thdot for p in self.pl]

        print(f"{(t/self.tf)*100:.1f}%")
        return thdot_data + thddot_data
    

    def integrate_system(self, t_int:np.ndarray) -> None:
        '''
        ### Integrates the system over the given time interval.
        Returns arrays consisting of each Pendulum object's:\n
        -angle data\n
        -angular vel data\n
        -angular accel data\n
        -x coor data\n
        -vel_x data\n
        -accel_x data\n
        -y coor data\n
        -vel_y data\n
        -accel_y data\n
        -kinetic energy data\n
        -potential energy data\n
        -total system energy data
        '''
        def _finite_central_difference(f):
            diff = np.zeros_like(f)
            diff[1:-1] = (f[2:] - f[:-2]) / (2*self.dt)
            return diff

        self.tf = t_int[-1]
        self.dt = t_int[1] - t_int[0]
        state0 = np.array([[p.th for p in self.pl], [p.thdot for p in self.pl]]).flatten()
        states = odeint(func=self._ode_system, y0=state0, t=t_int)

        x = lambda index: self.pl[index].L*np.sin(self.pl[index].th_data) + x(index-1) if index>=0 else 0
        y = lambda index: -self.pl[index].L*np.cos(self.pl[index].th_data) + y(index-1) if index>=0 else 0
        dxdt = lambda index: self.pl[index].L*self.pl[index].thdot_data*np.cos(self.pl[index].th_data) + dxdt(index-1) if index>=0 else 0
        dydt = lambda index: self.pl[index].L*self.pl[index].thdot_data*np.sin(self.pl[index].th_data) + dydt(index-1) if index>=0 else 0
        d2xdt2 = lambda index: self.pl[index].L * (self.pl[index].thddot_data*np.cos(self.pl[index].th_data) 
                                                   - self.pl[i].thdot_data**2*np.sin(self.pl[i].th_data)) + d2xdt2(index-1) if index>=0 else 0
        d2ydt2 = lambda index: self.pl[index].L * (self.pl[index].thddot_data*np.sin(self.pl[index].th_data) 
                                                   + self.pl[i].thdot_data**2*np.cos(self.pl[i].th_data)) + d2ydt2(index-1) if index>=0 else 0
        T = lambda index: 0.5*self.pl[index].m*(dxdt(index)**2 + dydt(index)**2)
        V = lambda index: self.pl[index].m*self.g*y(index)
        
        for i in range(self.n):
            self.pl[i].th_data = states.T[i]
            self.pl[i].thdot_data = states.T[i+self.n]
            self.pl[i].thddot_data = _finite_central_difference(self.pl[i].thdot_data)
            self.pl[i].x_data = x(i)
            self.pl[i].y_data = y(i)
            self.pl[i].xdot_data = dxdt(i)
            self.pl[i].ydot_data = dydt(i)
            self.pl[i].xddot_data = d2xdt2(i)
            self.pl[i].yddot_data = d2ydt2(i)
            self.pl[i].T = T(i)
            self.pl[i].V = V(i)
        
        self.E_tot = sum([self.pl[i].T + self.pl[i].V for i in range(self.n)])

        return [
            [self.pl[i].th_data for i in range(self.n)],
            [self.pl[i].thdot_data for i in range(self.n)],
            [self.pl[i].thddot_data for i in range(self.n)],
            [self.pl[i].x_data for i in range(self.n)],
            [self.pl[i].xdot_data for i in range(self.n)],
            [self.pl[i].xddot_data for i in range(self.n)],
            [self.pl[i].y_data for i in range(self.n)],
            [self.pl[i].ydot_data for i in range(self.n)],
            [self.pl[i].yddot_data for i in range(self.n)],
            [self.pl[i].T_data for i in range(self.n)],
            [self.pl[i].V_data for i in range(self.n)],
            self.E_tot
        ]
