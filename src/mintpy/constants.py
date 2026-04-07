############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, Feb 2022                           #
############################################################
# Recommend usage:
#   from mintpy.constants import SPEED_OF_LIGHT


SPEED_OF_LIGHT = 299792458  # meters per second

# Earth radius
# equatorial radius: a = 6378.1370e3
# polar      radius: b = 6356.7523e3
# arithmetic mean radius: R_1 = (2 * a + b) / 3 = 6371.0088e3
#   defined by IUGG and used in geophysics
EARTH_RADIUS = 6371.0088e3   # the arithmetic mean radius in meters


############################## Planetary Parameters ##############################
class PlanetaryBody():
    def __init__(self, name: str, radius: float, surface_gravity: float):
        self.name = name
        self.radius = radius                    # m
        self.surface_gravity = surface_gravity  # m/sec^2

Mercury = PlanetaryBody(name="Mercury", radius=2440e3,  surface_gravity=3.63)
Venus   = PlanetaryBody(name="Venus",   radius=6050e3,  surface_gravity=8.83)
Earth   = PlanetaryBody(name="Earth",   radius=6371e3,  surface_gravity=9.81)
Moon    = PlanetaryBody(name="Moon",    radius=1710e3,  surface_gravity=1.55)
Mars    = PlanetaryBody(name="Mars",    radius=3395e3,  surface_gravity=3.92)
Jupiter = PlanetaryBody(name="Jupiter", radius=71500e3, surface_gravity=25.9)
Saturn  = PlanetaryBody(name="Saturn",  radius=60000e3, surface_gravity=11.38)
