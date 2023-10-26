# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.

import matplotlib.pyplot as plt
import numpy as np
import math


PI                  = math.pi
KAPPA               = 6.674e-11     # N*m2/kg-2
MASS_EARTH          = 5.972e24      # kg
G                   = 9.81          # m/s-2
EARTH_RADIUS       = 6378e3         # meters


def calc_vector_grav_force(mass1, mass2, pos_vect):
    r_x = pos_vect[0]
    r_y = pos_vect[1]
    r_abs = np.sqrt(r_x**2 + r_y**2)
    grav_force = KAPPA * mass1 * mass2 * (-1) * np.array([r_x, r_y]) / r_abs**3
    return grav_force


class Vector:
    def __init__(self, x=0, y=0, z=0):
        self.x = x
        self.y = y
        self.z = z


class BodyObj:
    def __init__(self, mass, x_pos, y_pos, x_vel, y_vel):
        self.mass = mass
        self.position = np.array([x_pos, y_pos])
        self.velocity = np.array([x_vel, y_vel])

        self.force = np.array([0, 0])
        self.acceleration = self.force / self.mass

    def __calc_vector_force_from_other_body(self, other_body):
        r_x = self.position[0] - other_body.position[0]
        r_y = self.position[1] - other_body.position[1]
        r_abs = np.sqrt(r_x ** 2 + r_y ** 2)
        force_from_other_body = KAPPA * self.mass * other_body.mass * (-1) * np.array([r_x, r_y]) / r_abs ** 3
        return force_from_other_body

    def apply_forces_on_body(self, dt, other_bodies):
        force_from_other_bodies = 0
        for other_body in other_bodies:
            force_from_other_bodies += self.__calc_vector_force_from_other_body(other_body)

        self.force = force_from_other_bodies
        self.acceleration = self.force / self.mass
        dv_ball = self.acceleration * dt
        self.velocity = self.velocity + dv_ball
        ds_ball = self.velocity * dt
        self.position = self.position + ds_ball


def main():
    start_time = 0
    end_time = 100
    dt = 1e-3
    no_of_time_points = 1 + round((end_time - start_time) / dt)
    t_sim = np.linspace(start_time, end_time, no_of_time_points)
    dt = t_sim[1] - t_sim[0]
    x = []
    y = []

    inti_pos_obj_above_earth_surface = 25
    latitude_deg = 70
    lat_rad = 2 * PI / 360 * latitude_deg

    earth = BodyObj(MASS_EARTH, x_pos=0, y_pos=0, x_vel=0, y_vel=0)
    ball = BodyObj(mass=1, x_pos=(inti_pos_obj_above_earth_surface + EARTH_RADIUS) * np.sin(lat_rad),
                   y_pos=(inti_pos_obj_above_earth_surface + EARTH_RADIUS) * np.cos(lat_rad),
                   x_vel=0, y_vel=0)

    for t_tick in t_sim:
        earth.apply_forces_on_body(dt=dt, other_bodies=[ball])
        ball.apply_forces_on_body(dt=dt, other_bodies=[earth])

        x.append(ball.position[0])
        y.append(ball.position[1])
        if np.sqrt(ball.position[0] ** 2 + ball.position[1] ** 2) < EARTH_RADIUS:
            break

    for idx in range(len(y)):
        x[idx] = x[idx] - EARTH_RADIUS * np.sin(lat_rad)
        y[idx] = y[idx] - EARTH_RADIUS * np.cos(lat_rad)

    t = t_sim[:len(y)]

    print("t_end = " + str(t[-1]))
    print("axcel = " + str(ball.acceleration))
    # plt.plot(t, y)
    # plt.show()
    # plt.plot(t, x)
    # plt.show()
    plt.plot(x, y)
    plt.show()
    return


if __name__ == '__main__':
    print('PyCharm')
    main()
