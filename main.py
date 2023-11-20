# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.
import copy
import time

import matplotlib.pyplot as plt
import numpy as np
import math


PI                           = math.pi        # PI constant
KAPPA                        = 6.674e-11      # N*m2/kg-2

# data copied from https://nssdc.gsfc.nasa.gov/planetary/factsheet/
MERCURY_MASS                 = 0.330e24       # kg
MERCURY_PERIHELION_DISTANCE  = 46.0e9         # in meters distance from Sun in Mercury perihelion
MERCURY_PERIHELION_VELOTITY  = 58.97e3        # in m/s in Mercury perihelion

VENUS_MASS                   = 4.87e24        # kg
VENUS_PERIHELION_DISTANCE    = 107.5e9        # in meters distance from Sun in Venus perihelion
VENUS_PERIHELION_VELOTITY    = 35.02e3        # in m/s in Venus perihelion

EARTH_MASS                   = 5.97e24        # kg
EARTH_PERIHELION_DISTANCE    = 147.1e9        # in meters distance from Sun in Earth perihelion
EARTH_PERIHELION_VELOTITY    = 30.29e3        # in m/s in Earth perihelion
EARTH_SURFACE_RADIUS         = 6378e3         # meters

MARS_MASS                    = 0.642e24       # kg
MARS_PERIHELION_DISTANCE     = 206.7e9        # in meters distance from Sun in Mars perihelion
MARS_PERIHELION_VELOTITY     = 26.50e3        # in m/s in Mars perihelion

JUPITER_MASS                 = 1898e24        # kg
JUPITER_PERIHELION_DISTANCE  = 740.6e9        # in meters distance from Sun in Jupiter perihelion
JUPITER_PERIHELION_VELOTITY  = 13.72e3        # in m/s in Jupiter perihelion

SATURN_MASS                  = 568e24         # kg
SATURN_PERIHELION_DISTANCE   = 1357.6e9       # in meters distance from Sun in Saturn perihelion
SATURN_PERIHELION_VELOTITY   = 10.14e3        # in m/s in Saturn perihelion

URANUS_MASS                  = 86.8e24        # kg
URANUS_PERIHELION_DISTANCE   = 2732.7e9       # in meters distance from Sun in Uranus perihelion
URANUS_PERIHELION_VELOTITY   = 6.79e3         # in m/s in Uranus perihelion

NEPTUNE_MASS                 = 102e24         # kg
NEPTUNE_PERIHELION_DISTANCE  = 4471.1e9       # in meters distance from Sun in Neptune perihelion
NEPTUNE_PERIHELION_VELOTITY  = 5.47e3         # in m/s in Neptune perihelion

MOON_MASS                    = 0.073e24       # kg
MOON_PERIGEE_DISTANCE        = 0.3633e9       # in meters distance from Earth center in perigee
MOON_PERIGEE_VELOTITY        = 1.082e3        # m/s in Moon perigee from Earth POV

SUN_MASS                     = 1989000e24     # kg


MINUTE_sec             = 60
HOUR_sec               = 60 * MINUTE_sec
DAY_sec                = 24 * HOUR_sec
MONTH_sec              = 28 * DAY_sec
YEAR_sec               = 365 * DAY_sec


class Vector:
    def __init__(self, x=0, y=0):
        self.x = x
        self.y = y

    def __add__(self, other):
        if type(other) is type(self):
            x = self.x + other.x
            y = self.y + other.y
        else:
            x = self.x + other
            y = self.y + other
        return Vector(x, y)

    def __radd__(self, other):
        return self.__add__(other)

    def __mul__(self, other):
        if type(other) is type(self):
            x = self.x * other.x
            y = self.y * other.y
        else:
            x = self.x * other
            y = self.y * other
        return Vector(x, y)

    def __rmul__(self, other):
        return self.__mul__(other)

    def __sub__(self, other):
        return self.__add__(other * (-1))

    def __truediv__(self, other):
        return self.__mul__(1 / other)

    def __abs__(self):
        return np.sqrt(self.x ** 2 + self.y ** 2)


class BodyObj:
    def __init__(self, mass, x_pos, y_pos, x_vel, y_vel):
        self.mass = mass
        self.position = Vector(x_pos, y_pos)
        self.velocity = Vector(x_vel, y_vel)
        self.acceleration = Vector(0, 0)

    def precalc_vector_gravitational_acceleration_from_other_bodies(self, other_bodies):
        acceleration_from_other_bodies = Vector(0, 0)
        for other_body in other_bodies:
            if other_body is not self:
                r_vect = self.position - other_body.position
                acceleration_from_other_body = (-1) * KAPPA * other_body.mass * r_vect / abs(r_vect) ** 3
                acceleration_from_other_bodies += acceleration_from_other_body
        self.acceleration = acceleration_from_other_bodies

    def apply_acceleration_on_body(self, dt):
        dv = self.acceleration * dt
        self.velocity = self.velocity + dv
        ds = self.velocity * dt
        self.position = self.position + ds


def main():
    start_time = 0
    end_time = 3 * MONTH_sec
    dt = 1 * HOUR_sec
    no_of_time_points = 1 + round((end_time - start_time) / dt)
    t_sim = np.linspace(start_time, end_time, no_of_time_points)

    log_period = 1 * HOUR_sec
    stats_printout_period = min(0.01 * (end_time - start_time), 1 * YEAR_sec)

    sun = BodyObj(SUN_MASS,
                   x_pos=0, y_pos=0, x_vel=0, y_vel=0)
    mercury = BodyObj(mass=MERCURY_MASS,
                   x_pos=MERCURY_PERIHELION_DISTANCE, y_pos=0,
                   x_vel=0, y_vel=MERCURY_PERIHELION_VELOTITY)
    venus = BodyObj(mass=VENUS_MASS,
                   x_pos=VENUS_PERIHELION_DISTANCE, y_pos=0,
                   x_vel=0, y_vel=VENUS_PERIHELION_VELOTITY)
    earth = BodyObj(mass=EARTH_MASS,
                   x_pos=EARTH_PERIHELION_DISTANCE, y_pos=0,
                   x_vel=0, y_vel=EARTH_PERIHELION_VELOTITY)
    moon = BodyObj(mass=MOON_MASS,
                   x_pos=MOON_PERIGEE_DISTANCE + earth.position.x, y_pos=0 + earth.position.y,
                   x_vel=0+earth.velocity.x, y_vel=MOON_PERIGEE_VELOTITY + earth.velocity.y)
    mars = BodyObj(mass=MARS_MASS,
                   x_pos=MARS_PERIHELION_DISTANCE, y_pos=0,
                   x_vel=0, y_vel=MARS_PERIHELION_VELOTITY)
    jupiter = BodyObj(mass=JUPITER_MASS,
                   x_pos=JUPITER_PERIHELION_DISTANCE, y_pos=0,
                   x_vel=0, y_vel=JUPITER_PERIHELION_VELOTITY)
    saturn = BodyObj(mass=SATURN_MASS,
                   x_pos=SATURN_PERIHELION_DISTANCE, y_pos=0,
                   x_vel=0, y_vel=SATURN_PERIHELION_VELOTITY)
    uranus = BodyObj(mass=URANUS_MASS,
                   x_pos=URANUS_PERIHELION_DISTANCE, y_pos=0,
                   x_vel=0, y_vel=URANUS_PERIHELION_VELOTITY)
    neptune = BodyObj(mass=NEPTUNE_MASS,
                   x_pos=NEPTUNE_PERIHELION_DISTANCE, y_pos=0,
                   x_vel=0, y_vel=NEPTUNE_PERIHELION_VELOTITY)

    static_zero_point = BodyObj(mass=0,
                   x_pos=sun.position.x, y_pos=sun.position.y,
                   x_vel=0, y_vel=0)

    earth_balanced_with_moon = BodyObj(mass=EARTH_MASS,
                   x_pos=EARTH_PERIHELION_DISTANCE, y_pos=0,
                   x_vel=0, y_vel=-13.045*0)
    moon_balanced_with_earth = BodyObj(mass=MOON_MASS,
                   x_pos=MOON_PERIGEE_DISTANCE+earth_balanced_with_moon.position.x, y_pos=0+earth_balanced_with_moon.position.y,
                   x_vel=0, y_vel=MOON_PERIGEE_VELOTITY+earth_balanced_with_moon.velocity.y)

    static_zero_point = BodyObj(mass=0,
                   x_pos=earth_balanced_with_moon.position.x, y_pos=earth_balanced_with_moon.position.y,
                   x_vel=0, y_vel=0)

    # TODO create a system class. At construction, it must verify if two objects overlap
    # solar_system_bodies = [sun, earth, jupiter, saturn]
    # solar_system_bodies = [sun, venus, earth, moon]
    # solar_system_bodies = [sun, mercury, venus, earth, mars, jupiter, saturn, uranus, neptune, moon]
    # solar_system_bodies = [sun, earth]
    solar_system_bodies = [earth_balanced_with_moon, moon_balanced_with_earth]
    plotted_body = earth_balanced_with_moon
    observer_body = static_zero_point

    log_plotted_body_data = dict()
    log_observer_data = dict()

    t1 = time.time()

    for t_tick in t_sim:
        for body in solar_system_bodies:
            body.precalc_vector_gravitational_acceleration_from_other_bodies(other_bodies=solar_system_bodies)
        for body in solar_system_bodies:
            body.apply_acceleration_on_body(dt=dt)

        if t_tick % log_period == 0:
            log_plotted_body_data.update({t_tick: copy.copy(plotted_body)})
            log_observer_data.update({t_tick: copy.copy(observer_body)})

        if t_tick % stats_printout_period == 0 and t_tick > 0:
            tx = time.time()
            no_of_simul_days = t_tick / DAY_sec
            no_of_simul_years = t_tick / YEAR_sec
            complete_percent = t_tick / end_time * 100

            if stats_printout_period > 100 * DAY_sec:
                printout_message = \
                    "{:5.2f}% of simulation completed. {: 7d} years simulated by speed {:5.2f} years per second"\
                        .format(complete_percent,
                                round(no_of_simul_years),
                                no_of_simul_years / (tx - t1))
            else:
                printout_message = \
                    "{:5.2f}% of simulation completed. {: 7d} days simulated by speed {:5.0f} days per second" \
                        .format(complete_percent,
                                round(no_of_simul_days),
                                no_of_simul_days / (tx - t1))

            progress_bar = "  ||" + round(complete_percent) * "=" + (100 - round(complete_percent)) * " " + "||"
            printout_message += progress_bar

            clear_message = "".join(['\b'] * len(printout_message))
            print("{:s}{:s}".format(clear_message, printout_message), end='', flush=True)

    t2 = time.time()
    print("\nSimulation took " + str(round(t2 - t1, 2)) + " seconds")
    # ==============================================================
    # ==============================================================
    print('\nData post-processing')
    t1 = time.time()
    analyzed_distance = []
    for body, observer in zip(log_plotted_body_data.values(), log_observer_data.values()):
        dist_from_observer = abs(body.position - observer.position)
        analyzed_distance.append(round(dist_from_observer))

    simout_planet_aphelion = max(analyzed_distance)
    simout_planet_perihelion = min(analyzed_distance)
    print("Planet aphelion: " + str(round(simout_planet_aphelion / 1e9, 3)) + " mil km")
    print("Planet perihelion: " + str(round(simout_planet_perihelion / 1e9, 3)) + " mil km")

    log_to_file = []
    for t_stamp, body in log_plotted_body_data.items():
        log_line = [str(item) for item in (t_stamp, body.position.x, body.position.y, '\n')]
        log_to_file.append(str.join(',', log_line))
    with open("logfile.txt", "w") as fw:
        fw.writelines(log_to_file)

    t2 = time.time()
    print("\nData post-processing took " + str(round(t2 - t1, 2)) + " seconds")
    # ==============================================================
    # ==============================================================
    print('\nData imaging')
    track = [body.position - obsv.position for body, obsv in zip(log_plotted_body_data.values(), log_observer_data.values())]
    x = [item.x for item in track]
    y = [item.y for item in track]
    plt.plot(x, y)
    plt.show()
    return


if __name__ == '__main__':
    print('N-body simulation of solar system')
    main()
