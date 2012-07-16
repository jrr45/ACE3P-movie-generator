#!/usr/bin/env python

from paraview.simple import *
from ace3p_animation import Ace3pAnimataion
from general_camera import GeneralCamera, CameraMovement
from math import pi
import glob

################################################################################
# set up world
anime = Ace3pAnimataion(mesh_file = '/home/jrr45/gun_current/gun.ncdf',
            mode_files=['/home/jrr45/gun_current/omega3p_results/omega3p.l0.m0000.1.4950919e+09.mod'], 
            trac3p_results="",
            reflections=['X'],
            sidesetIDs=[4,6],
            field='efield')

################################################################################
# set up camera
x, y, z = anime.get_center()
time_steps=anime.get_timestep_values()

camera = GeneralCamera(camera_movements=[],
                       focus=[x, y, z], 
                       radius=.83, theta=-pi/2, phi=43*pi/100,
                       up=[-1.0,0.0,0.0], psi=0)

################################################################################
# run animation
anime.play(camera=camera, 
           background_color=[0.0, 0.0, 0.0],
           time_steps=time_steps, antialias=2)