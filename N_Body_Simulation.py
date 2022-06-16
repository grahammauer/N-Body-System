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
    forcevectors = []

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
        forcevectors.append(vector)

    # return the list of force vectors out of the function
    return forcevectors
