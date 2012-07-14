#!/usr/bin/env python

################################################################################

try: paraview.simple
except: from paraview.simple import *
from ace3p_animation import Ace3pAnimataion
from general_camera import GeneralCamera, CameraMovement
from math import pi
import glob

################################################################################
root = "/home/jrr45/ptesting/taper/"
anime = Ace3pAnimataion(mesh_file = '/home/jrr45/gun_current/gun.ncdf',#mesh_file = (root + "taper.ncdf"),
            mode_files = ['/home/jrr45/gun_current/omega3p_results/omega3p.l0.m0000.1.4950919e+09.mod'], #mode_files = glob.glob(root+"/t3p_results/OUTPUT/mymonts_t[0-9]*ps.out.mod"),#trac3p_results = "/home/jrr45/small_gap/track3p_results_6.0/",
            reflections=['X'],
            sidesetIDs=[4,6],
            field='efield')

################################################################################
#find the center of the mesh
x, y, z = anime.get_center()

camera = GeneralCamera(camera_movements=[],
                       focus=[x, y, z], 
                       radius=.83, theta=-pi/2, phi=43*pi/100,
                       up=[-1.0,0.0,0.0], psi=0)

anime.play(camera=camera, background_color=[0.0, 0.0, 0.0],
           GUI=False, antialias=2)