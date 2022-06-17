#######################################################################
# Modules
#######################################################################

import numpy as np

#######################################################################
# Physics
#######################################################################
def forceMagnitude(mi, mj, sep):
    """
    Compute magnitude of gravitational force between two particles.
    """
    G = 6.67e-11                # m3 kg-1 s-2
    return G * mi * mj / sep**2 # N

def magnitude(vec):
    """
    Compute magnitude of any vector with an arbitrary number of elements.
    """
    return np.sqrt(np.sum(vec**2))

def unitDirectionVector(pos_a, pos_b):
    """
    Create unit direction vector from pos_a to pos_b
    """

    # calculate the separation between the two vectors
    separation = pos_b - pos_a

    # divide vector components by vector magnitude to make unit vector
    return separation/magnitude(separation)

def forceVector(mi, mj, pos_i, pos_j):
    """
    Compute gravitational force vector exerted on particle i by particle j.
    """

    # compute the magnitude of the distance between positions
    distance = magnitude(pos_i - pos_j)

    # compute the magnitude of the force
    force = forceMagnitude(mi, mj, distance)

    # calculate the unit direction vector of the force
    direction = unitDirectionVector(pos_i, pos_j)

    return force*direction # a numpy array, with units of Newtons

def calculateForceVectors(masses, positions):
    """
    Compute net gravitational force vectors on particles,
    given a list of masses and positions for all of them.
    """

    # how many particles are there?
    N = len(positions)

    # create an empty list, which we will fill with force vectors
    forcevectors = np.zeros_like(positions)

    # loop over particles for which we want the force vector
    for i in range(N):

        # create a force vector with all three elements as zero
        vector = np.zeros(3)

        # loop over all the particles we need to include in the force sum
        for j in range(N):

            # as long as i and j are not the same...
            if j != i:

                # ...add in the force vector of particle j acting on particle i
                vector += forceVector(masses[i], masses[j], positions[i], positions[j])

        # append this force vector into the list of force vectors
        forcevectors[i] = vector

    # return the list of force vectors out of the function
    return forcevectors

def calculateAccelerationVectors(positions, masses):
    """
    Compute net acceleration vectors on particles,
    given a list of masses and net accelerations for all of them.
    """
    forcevectors = calculateForceVectors(masses, positions)

    accelerations = np.zeros_like(forcevectors)
    for i in range(len(forcevectors)):
        accelerations[i] = forcevectors[i] / masses[i]

    return accelerations

#######################################################################
# Propogation Through Time
#######################################################################

def updateParticles(masses, positions, velocities, dt):
    '''
    Evolve particles in time via leap-frog integrator scheme. This function
    takes masses, positions, velocities, and a time step dt as
    
    Parameters
    ----------
    masses : np.ndarray
        1-D array containing masses for all particles, in kg
        It has length N, where N is the number of particles.
    positions : np.ndarray
        2-D array containing (x, y, z) positions for all particles.
        Shape is (N, 3) where N is the number of particles.
    velocities : np.ndarray
        2-D array containing (x, y, z) velocities for all particles.
        Shape is (N, 3) where N is the number of particles.
    dt : float
        Evolve system for time dt (in seconds).
    
    Returns
    -------
    Updated particle positions and particle velocities, each being a 2-D
    array with shape (N, 3), where N is the number of particles.
    '''

    positionsNew = positions + velocities * dt

    velocitiesNew = velocities + calculateAccelerationVectors(positions, masses) * dt

    return positionsNew, velocitiesNew

def calculateTrajectories(masses, positions, velocities, steps, dt):
    '''
    Parameters
    ----------
    masses : np.ndarray
        a 1D array with nParticle elements that has the masses for each 
        particle
    positions : np.ndarray
        a 2D array of shape (nParticles, nDimensions) that contains a 
        (x,y,z) intitial position for each particle
    velocities : np.ndarray
        a 2D array of shape (nParticles, nDimensions) that contains a 
        (x,y,z) initial velocity for each particle
    steps : integer
        An integer number of dt steps for the function to go through
    dt : float
        Time spacing between steps
        
    
    Returns
    -------
    times : np.ndarray
        a 1D array of all the times that the function calculated positions
        and velocities for
    positions : np.ndarray
        a 3D array of shape (times, nParticles, nDimensions) containing the
        positions of each particle at every time.
    velocities : np.ndarray
        a 3D array of shape (times, nParticles, nDimensions) containing the
        velocities of each particle at every time.
    '''
    # Create the array of times
    times = np.arange(0, steps+1) * dt
    
    # Create blank position and velocity arrays to iterate through
    positionsArray = np.array([positions])
    velocitiesArray = np.array([velocities])
    
    # Iterate through each time step with time dt
    for i in range(len(times)):
        
        # Run updateParticles()
        leapfrog = updateParticles(masses,positions,velocities,dt)
        
        # Save the outputs
        positionsArray = np.append(positionsArray, [leapfrog[0]])
        velocitiesArray = np.append(velocitiesArray, [leapfrog[1]])       
        
        # Save outputs to put back into updateParticles()
        positions = leapfrog[0]
        velocities = leapfrog[1]
    
    # Return final arrays, I needed to reshape the positions and velocities
    # into the correct format because of how I wrote the above code that appends.
    return positionsArray.reshape(steps+2,len(masses),3), velocitiesArray.reshape(steps+2,len(masses),3), times
