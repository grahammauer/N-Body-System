# N-Body-System

Simulates a gravitational N-Body system in Python3. Objects are inputted as a collection of numpy arrays. To define an object; initial mass, position, and velocity must be known. From there, two or more object interact through a leapfrog scheme where they advance through time by small steps and update their position and velocity. 

The simulation is ran by the `calculateTrajectories()` function.
