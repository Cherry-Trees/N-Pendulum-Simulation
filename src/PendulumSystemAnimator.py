import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.animation import FuncAnimation
import numpy as np


class PendulumSystemAnimator:

    BOB_SIZE_SCALAR = 26.66
    AX_BOUNDS_SCALAR = 1.2
    TRAIL_DENSITY = 20
    TRAIL_LENGTH = 2
    TRAIL_COLOR_SENSITIVITY = 2.33

    @staticmethod
    def animate(pendulumSystemList:list, t, dt, figsize=(6.4, 6.4), bgcolor="white"):
        rcParams["toolbar"] = "None"
        DEFAULT_COLORS = rcParams["axes.prop_cycle"].by_key()["color"]
        fig = plt.figure(figsize=figsize, facecolor=bgcolor)
        ax = fig.add_axes((0, 0, 1, 1), facecolor=bgcolor)

        BOUND = PendulumSystemAnimator.AX_BOUNDS_SCALAR*max([ps.L_tot for ps in pendulumSystemList])
        ax.set_xlim(-BOUND, BOUND)
        ax.set_ylim(-BOUND, BOUND)
        t_int = np.arange(0, t, dt)

        for ps in pendulumSystemList:
            ps.integrate_system(t_int)
            for i in range(ps.n):
                if ps.pl[i].visualstyle == "trail":
                    ps.pl[i].trail = [ax.plot([], [], lw=1, c="white", alpha=0, solid_capstyle="butt")[0]
                                for _ in range(PendulumSystemAnimator.TRAIL_DENSITY)]
                
                if ps.pl[i].visualstyle == "tracer":
                    ps.pl[i].tracer, = ax.plot([],[], color=ps.pl[i].visualstylecolors[0] if ps.pl[i].visualstylecolors != "DEFAULT" else DEFAULT_COLORS[i%len(DEFAULT_COLORS)])
                
                if ps.pl[i].visualstyle == "vector":
                    ps.pl[i].v_vector = ax.annotate("", xy=(0,0), xytext=(0, 0), arrowprops={"facecolor": "lightblue"})
                    ps.pl[i].a_vector = ax.annotate("", xy=(0,0), xytext=(0, 0), arrowprops={"facecolor": "coral"})

            for i in range(ps.n):
                ps.pl[i].rod_sprite, = ax.plot([],[], linewidth=ps.pl[i].rodsize, 
                                               color=ps.pl[i].rodcolor if ps.pl[i].rodcolor != "DEFAULT" else DEFAULT_COLORS[i%len(DEFAULT_COLORS)])
                if ps.pl[i].showbob:
                    ps.pl[i].bob_sprite, = ax.plot([],[], "o", markersize=np.log(ps.pl[i].m+2)*(PendulumSystemAnimator.BOB_SIZE_SCALAR/np.sqrt(2*BOUND)), 
                                                   color=ps.pl[i].bobcolor if ps.pl[i].bobcolor != "DEFAULT" else DEFAULT_COLORS[i%len(DEFAULT_COLORS)])
        
        def _update(frame) -> list:
            obj_array = []
          
            for ps in pendulumSystemList:
                for i in range(ps.n):
                    if ps.pl[i].visualstyle == "trail":
                        trail_colors = ps.pl[i].visualstylecolors if ps.pl[i].visualstylecolors != "DEFAULT" else [DEFAULT_COLORS[i%len(DEFAULT_COLORS)]]
                        for segment in range(PendulumSystemAnimator.TRAIL_DENSITY):
                            obj_array.append(ps.pl[i].trail[segment])
                            frame_min = max(frame - (PendulumSystemAnimator.TRAIL_DENSITY-segment)*PendulumSystemAnimator.TRAIL_LENGTH, 0)
                            frame_max = frame_min + PendulumSystemAnimator.TRAIL_LENGTH + 1
                            segment_alpha = (segment/PendulumSystemAnimator.TRAIL_DENSITY)**3

                            ps.pl[i].trail[segment].set_data(ps.pl[i].x_data[frame_min:frame_max], ps.pl[i].y_data[frame_min:frame_max])
                            ps.pl[i].trail[segment].set_alpha(segment_alpha)
                            ps.pl[i].trail[segment].set_color(trail_colors[int(min(abs(ps.pl[i].thdot_data[frame_min])/PendulumSystemAnimator.TRAIL_COLOR_SENSITIVITY,
                                                                            len(trail_colors)-1))])
                    
                    if ps.pl[i].visualstyle == "tracer":
                        obj_array.append(ps.pl[i].tracer)
                        ps.pl[i].tracer.set_data(ps.pl[i].x_data[:frame], ps.pl[i].y_data[:frame])
                    
                    if ps.pl[i].visualstyle == "vector":
                        obj_array.append(ps.pl[i].v_vector)
                        obj_array.append(ps.pl[i].a_vector)
                        ps.pl[i].v_vector.xy = ps.pl[i].x_data[frame] + ps.pl[i].xdot_data[frame]/6, ps.pl[i].y_data[frame] + ps.pl[i].ydot_data[frame]/6
                        ps.pl[i].a_vector.xy = ps.pl[i].x_data[frame] + ps.pl[i].xddot_data[frame]/36, ps.pl[i].y_data[frame] + ps.pl[i].yddot_data[frame]/36
                        
                        ps.pl[i].v_vector.set_position([ps.pl[i].x_data[frame], ps.pl[i].y_data[frame]])
                        ps.pl[i].a_vector.set_position([ps.pl[i].x_data[frame], ps.pl[i].y_data[frame]])
                
                for i in range(ps.n):
                    ps.pl[i].rod_sprite.set_data([(ps.pl[i-1].x_data[frame:frame+1]) if i-1>-1 else [0], ps.pl[i].x_data[frame:frame+1]],
                                            [ps.pl[i-1].y_data[frame:frame+1] if i-1>-1 else [0], ps.pl[i].y_data[frame:frame+1]])
                    obj_array.append(ps.pl[i].rod_sprite)
                for i in range(ps.n):
                    if ps.pl[i].showbob:
                        ps.pl[i].bob_sprite.set_data(ps.pl[i].x_data[frame:frame+1], ps.pl[i].y_data[frame:frame+1])
                        obj_array.append(ps.pl[i].bob_sprite)

            return obj_array
        
        ANIMATION = FuncAnimation(fig, _update, len(t_int), interval=1, blit=True)
        plt.show()
