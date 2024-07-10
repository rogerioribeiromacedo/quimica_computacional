# -*- coding: utf-8 -*-
"""
Newton's Law.

@author: rogerio
"""


class Particle:
    """Class Particle"""

    def __init__(self, mass, velocity=0):
        self.mass = mass
        self.velocity = velocity

    def apply_a_force(self, force, time):
        """
        Apply a force over the particle.

        Parameters
        ----------
        force : float
            Force in Newton.
        time : int
            Time in second.

        Returns
        -------
        None.

        """
        # Calculate the acceleration of the particle
        self.acceleration = force / self.mass

        # Update de velocity of the particle
        self.velocity += self.acceleration * time


if __name__ == "__main__":
    # Create a particle with mass equal to 5 kg
    parti = Particle(5)
    print(f"Initial velocity: {parti.velocity}")

    # Apply a force of 10 N for 1 second.
    parti.apply_a_force(20, 1)
    print(f"Velocity: {parti.velocity}")
