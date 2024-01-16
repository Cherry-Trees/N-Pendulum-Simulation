# Description

An n-pendulum simulation written in Python 3 that utilizes Scipy's odeint function for integration and Matplotlib for visualization. The program solves a set of n second-order differential equations derived from the Lagrangian formulation of an n-pendulum system.

# Examples

### Standard Double Pendulum
A double pendulum (n=2) with both masses 1, rod lengths 1, and initial angles 90 deg yields the following results:
<p align="center">
  

  <img src="https://github.com/Cherry-Trees/n-pendulum/blob/main/examples/pendulum2.gif" width="400" height="400" />
  <img src="https://github.com/Cherry-Trees/n-pendulum/blob/main/examples/pendulum2vectors.gif" width="400" height="400" />
    


</p>


Additionally, we can plot some of the system's properties, including each of the pendulum's angles over time, the path the pendulums took (x and y), the pendulum's individual energies, and the total energy drift of the system:
<p align="center">
  <img src="https://github.com/Cherry-Trees/n-pendulum/blob/main/examples/angle2.png" width="400" height="300"/>
  <img src="https://github.com/Cherry-Trees/n-pendulum/blob/main/examples/path2.png" width="400" height="300"/>
  <img src="https://github.com/Cherry-Trees/n-pendulum/blob/main/examples/energy2.png" width="400" height="300"/>
  <img src="https://github.com/Cherry-Trees/n-pendulum/blob/main/examples/drift2.png" width="400" height="300"/>
  
</p>













