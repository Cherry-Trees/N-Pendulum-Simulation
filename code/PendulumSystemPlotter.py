import matplotlib.pyplot as plt
import numpy as np


class PendulumSystemPlotter:

    @staticmethod
    def plot(pendulumSystem, t, dt, figsize=(10, 7.5)):
        fig = plt.figure(figsize=figsize)
        ax = fig.add_axes((0.1, 0.1, 0.75, 0.8))
        t_int = np.arange(0, t, dt)
        pendulumSystem.integrate_system(t_int)
        for i in range(pendulumSystem.n):
            ax.plot(t_int, pendulumSystem.pl[i].th_data, label=rf"$P_{i+1}$")

        ax.set_title(f"Angle v. Time")
        ax.set_xlabel("t")
        ax.set_ylabel(r"${\theta}$(t)")
        ax.legend(bbox_to_anchor=(1.135, 1), ncol=1)

        global plotNum
        plotNum = 0
        
        def _on_click(event):
            global plotNum
            plotNum = (plotNum+1) % 4

            event.canvas.figure.clear()
            ax = fig.add_axes((0.1, 0.1, 0.75, 0.8))
            ax.set_title(f"n={pendulumSystem.n}")
  
            if plotNum==0:
                for i in range(pendulumSystem.n):
                    ax.plot(t_int, pendulumSystem.pl[i].th_data, label=rf"$P_{i+1}$")
                ax.set_title("Angle v. Time")
                ax.set_xlabel("t")
                ax.set_ylabel(r"${\theta}$(t)")
                ax.legend(bbox_to_anchor=(1.135, 1), ncol=1)
                ax.figure.canvas.draw()
            
            if plotNum==1:
                ax.set_xlim(-pendulumSystem.L_tot*1.2*(7.5/5.5), pendulumSystem.L_tot*1.2*(7.5/5.5))
                ax.set_ylim(-pendulumSystem.L_tot*1.2, pendulumSystem.L_tot*1.2)

                for i in range(pendulumSystem.n):
                    ax.plot(pendulumSystem.pl[i].x_data, pendulumSystem.pl[i].y_data, label=rf"$P_{i+1}$")
                ax.set_title("Position")
                ax.set_xlabel("x(t)")
                ax.set_ylabel("y(t)")
                ax.legend(bbox_to_anchor=(1.135, 1), ncol=1)
                ax.figure.canvas.draw()
            
            if plotNum==2:
                for i in range(pendulumSystem.n):
                    ax.plot(t_int, pendulumSystem.pl[i].T + pendulumSystem.pl[i].V, label=rf"$P_{i+1}$")
                ax.set_title("Energy v. Time")
                ax.set_xlabel("t")
                ax.set_ylabel("E(t)")
                ax.legend(bbox_to_anchor=(1.135, 1), ncol=1)
                ax.figure.canvas.draw()
            
            if plotNum==3:
                ax.plot(t_int, pendulumSystem.E_tot, label=r"$E_{tot}$")
                ax.set_title("Total Energy Drift")
                ax.set_xlabel("t")
                ax.set_ylabel("E(t)")
                ax.legend(bbox_to_anchor=(1.135, 1), ncol=1)
                ax.figure.canvas.draw()

        fig.canvas.mpl_connect('button_press_event', _on_click)
        plt.show()
