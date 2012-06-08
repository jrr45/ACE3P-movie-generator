#!/usr/bin/env python
##!/usr/bin/pvbatch --use-offscreen-rendering

################################################################################

try: paraview.simple
except: from paraview.simple import *
from ace3p_animation import Ace3pAnimataion
from general_camera import GeneralCamera, CameraMovement
from math import pi

################################################################################
anime = Ace3pAnimataion(mesh_file = '/home/jrr45/small_gap/small_gap_6.0.ncdf',
            mode_file = '/home/jrr45/small_gap/omega3p_results_6.0/omega3p.l0.m0000.1.4957339e+09.mod',
            trac3p_results = "/home/jrr45/small_gap/track3p_results_6.0/",
            reflections=['X'])

################################################################################
#find the center of the mesh
x, y, z = anime.get_center()

camera = GeneralCamera(camera_movements=[],
                          focus=[x, y, z], 
                          radius=.83, theta=-pi/2, phi=43*pi/100,
                          up=[-1.0,0.0,0.0], psi=0)

anime.play(camera=camera, 
           GUI=False, antialias=2)