# -*- coding: utf-8 -*-
"""
Simulating particles that interact with Lennard-Jones potential.

That is a simple example of how we can simulate interactions between
particles. It does not pretend to be a complete molecular dynamics program,
so we make some assumptions.
  - The particles are represented by circles with the same radius.
  - Every particle will start at a random position and with a random velocity.
  - The particles interact with each other through a Lennard-Jones potential.
  - Use of periodic boundary conditions.

Lennard-Jones equaion:

    U (r) = (4 x epsilon) x [ (sigma/r)**12 - (sigma/r)**6 ]

--------------------------
"Kiti" is a word derived from Nheengatu and functions as a postposition,
meaning "to" or "towards".
    - Asu makiti aputári.
--------------------------

Author...........: Rogério Ribeiro Macêdo.
Curriculum Lattes: http://lattes.cnpq.br/8806221981552346
Last update......: July 12th, 2024.
"""
import sys
import numpy as np

# Length of text string
n_ljust = 50

# Radius for particle
radius_particle = 0.3


def head_msg():
    """
    Header message.

    Returns
    -------
    None.

    """
    print()
    print("-".center(79, "-"))
    print(f'{"|":<1} {"UNIFEI - Federal University of Itajubá":^75} '
          f'{"|":>1}')
    print(f'{"|":<1} {"LaQC - Computational Chemistry Laboratory":^75} '
          f'{"|":>1}')
    print(f'{"|":<1} {" ":^75} {"|":>1}')
    print(f'{"|":<1} {">>> use the command [exit] to terminate the program <<<":^75} '
          f'{"|":>1}')
    print(f'{"|":<1} {" ":^75} {"|":>1}')
    print("-".center(79, "-"))
    print(f'{"|":<1} {"Simulating particles that interact with Lennard-Jones potential.":^75} '
          f'{"|":>1}')
    print("-".center(79, "-"))
    print()


def tchau():
    """
    We are leaving and saying goodbye (tchau, inté!!).

    Returns
    -------
    None.

    """
    print("")
    print("-".center(79, "-"))
    print(f'{"|":<1} {"UNIFEI - Federal University of Itajubá":^75} '
          f'{"|":>1}')
    print(f'{"|":<1} {"LaQC - Computational Chemistry Laboratory":^75} '
          f'{"|":>1}')
    print(f'{"|":<1} {" ":^75} {"|":>1}')
    print(f'{"|":<1} {"Inté!!!!":^75} '
          f'{"|":>1}')
    print(f'{"|":<1} {" ":^75} {"|":>1}')
    print("-".center(79, "-"))
    print("")
    sys.exit()


def valid_number(val, default, value_type="integer"):
    """
    Valid the number.

    Parameters
    ----------
    val : str
        Value.
    default : int/float
        Value default.
    value_type : str, optional
        Type of number. The default is "integer".

    Returns
    -------
    value : int
        If the number is valid or not.

    """
    try:
        if len(val) == 0:
            val = default
        else:
            if value_type == "integer":
                val = int(val)
            elif value_type == "float":
                val = float(val)
    except ValueError:
        print(f" - The value must be an(a) {value_type}.")
        val = 0

    return val


def read_number_particles():
    try:
        number_particles = input("Number of particless "
                                 "[100]".ljust(n_ljust, ".") + ": ").strip()
        if number_particles == "exit":
            tchau()
        else:
            if len(number_particles) > 0:
                number_particles = int(number_particles)
            else:
                # Default value
                number_particles = 100
    except ValueError:
        print(" - The value must be an integer value.")
        return 0

    return number_particles


def read_lenght_box():
    try:
        lenght_box = input("Lenght of box "
                           "[10.0]".ljust(n_ljust, ".") + ": ").strip()
        if lenght_box == "exit":
            tchau()
        else:
            if len(lenght_box) > 0:
                lenght_box = float(lenght_box)
            else:
                # Default value
                lenght_box = 10.0
    except ValueError:
        print(" - The value must be an float value.")
        return 0

    return lenght_box


def read_duration_simul():
    try:
        duration_simul = input("Duration of simulation "
                               "[10]".ljust(n_ljust, ".") + ": ").strip()
        if duration_simul == "exit":
            tchau()
        else:
            if len(duration_simul) > 0:
                duration_simul = int(duration_simul)
            else:
                # Default value
                duration_simul = 10
    except ValueError:
        print(" - The value must be an integer value.")
        return 0

    return duration_simul


def read_number_steps():
    try:
        number_steps = input("Number of steps "
                             "[10]".ljust(n_ljust, ".") + ": ").strip()
        if number_steps == "exit":
            tchau()
        else:
            if len(number_steps) > 0:
                number_steps = int(number_steps)
            else:
                # Default value
                number_steps = 10
    except ValueError:
        print(" - The value must be an integer value.")
        return 0

    return number_steps


def read_initial_velocity():
    try:
        initial_velocity = input("Initial velocity "
                                 "[1.5]".ljust(n_ljust, ".") + ": ").strip()
        if initial_velocity == "exit":
            tchau()
        else:
            if len(initial_velocity) > 0:
                initial_velocity = float(initial_velocity)
            else:
                # Default value
                initial_velocity = 1.5
    except ValueError:
        print(" - The value must be an float value.")
        return 0

    return initial_velocity


def make_initial_structure(number_particles, lenght_box):
    """
    Construct a gro file with the initial structure.

    Parameters
    ----------
    number_particles : int
        Number of particle in a box.
    lenght_box : float
        Lenght of a box.

    Returns
    -------
    None.

    """
    # Create enough positions in grid for the particles
    grid_size = int(np.ceil(np.sqrt(number_particles)))
    print(grid_size)
    with open("initial_structure.gro", "w") as f_gro:
        f_gro.write("Simulating particles that interact with Lennard-Jones potential.\n")
        f_gro.write("{:>5d}\n".format(number_particles))
        f_gro.write("{:>5d}{:<5}{:>5}{:>5}{:>8.3f}{:>8.3f}{:>8.3f}\n".format(1, "PAR", "H", 1, 1.260, 1.251, 1.290))
        f_gro.write("{:>10.5f}{:>10.5f}{:>10.5f}\n".format(lenght_box, lenght_box, lenght_box))

    f_gro.close()


def simulation(number_particles, lenght_box, duration_simul, number_steps, initial_velocity):
    """
    Make simulation.

    Parameters
    ----------
    number_particles : TYPE
        DESCRIPTION.
    lenght_box : TYPE
        DESCRIPTION.
    duration_simul : TYPE
        DESCRIPTION.
    number_steps : TYPE
        DESCRIPTION.
    initial_velocity : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    pass


def main():
    """
    Principal function.

    Returns
    -------
    None.

    """
    number_particles = read_number_particles()
    if number_particles:
        lenght_box = read_lenght_box()
        if lenght_box:
            duration_simul = read_duration_simul()
            if duration_simul:
                number_steps = read_number_steps()
                if number_steps:
                    initial_velocity = read_initial_velocity()
                    if initial_velocity:
                        # Initial structure
                        make_initial_structure(int(number_particles), lenght_box)

                        # Simulation
                        simulation(int(number_particles), int(lenght_box),
                                   int(duration_simul), int(number_steps),
                                   float(initial_velocity))

                        # print(number_particles, lenght_box, duration_simul, number_steps, initial_velocity)


if __name__ == "__main__":
    # Show header message.
    head_msg()

    # Main
    main()
