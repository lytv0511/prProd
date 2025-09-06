**Woodberry Forest School's 2026 USAYPT Problem 4 Simulation**
===========================================================

**Multi-Bounce Kinematics Simulation**
------------------------------------

This repository contains a Python-based simulation for the 2026 USAYPT Problem 4 on Multi-Bounce Kinematics. The simulation models the motion of a ball bouncing on a surface, taking into account aerodynamic forces, gravity, and collisions with the surface.

**Table of Contents**
-----------------

1. [Introduction](#introduction)
2. [Simulation Overview](#simulation-overview)
3. [Code Structure](#code-structure)
4. [Usage](#usage)
5. [Contributing](#contributing)
6. [License](#license)

**Introduction**
---------------

The 2026 USAYPT Problem 4 on Multi-Bounce Kinematics challenges students to model the motion of a ball bouncing on a surface, considering various factors such as aerodynamic forces, gravity, and collisions. This repository provides a Python-based solution to this problem, utilizing numerical methods to simulate the ball's motion.

**Simulation Overview**
---------------------

The simulation consists of the following components:

* **Ball**: Represents the ball's physical properties, such as mass, radius, and moment of inertia.
* **Aerodynamics**: Models the aerodynamic forces acting on the ball, including drag and lift.
* **Surface**: Represents the surface on which the ball bounces, including its normal vector and size.
* **Collision Detection**: Detects collisions between the ball and the surface, updating the ball's velocity and angular velocity accordingly.
* **Numerical Integration**: Utilizes the Runge-Kutta method to numerically integrate the ball's motion over time.

**Code Structure**
-----------------

The code is organized into the following modules:

* `forces.py`: Contains functions for calculating aerodynamic forces and torques.
* `bouncesim.py`: Defines the `Simulator` class, which manages the simulation and updates the ball's state.
* `main.py`: Provides a simple example usage of the simulation.

**Usage**
---------

To run the simulation, simply execute the `main.py` script. This will simulate the ball's motion and plot its trajectory.

**Contributing**
--------------

Contributions are welcome! If you'd like to improve the simulation or add new features, please fork the repository and submit a pull request.

**License**
-------

This repository is licensed under the MIT License. See `LICENSE` for details.