import numpy as np
from PendulumSystem import PendulumSystem as PS
from PendulumSystemAnimator import PendulumSystemAnimator as PSA
from PendulumSystemPlotter import PendulumSystemPlotter as PSP

ps = PS()
TRAIL_COLORS = [
            "firebrick",
            "red",
            "orangered",
            "orange",
            "gold",
            "yellow",
            "lightyellow",
            "white"
         ]

ps.addPendulum(m=1, L=1, th0=np.radians(90), bobcolor="white", rodcolor="white", rodsize=0.5, visualstylecolors=TRAIL_COLORS)
ps.addPendulum(m=1, L=0.75, th0=np.radians(90), bobcolor="white", rodcolor="white", rodsize=0.5, visualstylecolors=TRAIL_COLORS)
ps.addPendulum(m=1, L=0.5, th0=np.radians(120), bobcolor="white", rodcolor="white", rodsize=0.5, visualstylecolors=TRAIL_COLORS)

PSA.animate(pendulumSystemList=[ps], t=20, dt=0.025, bgcolor="black")
# PSP.plot(pendulumSystem=ps, t=20, dt=0.025)
