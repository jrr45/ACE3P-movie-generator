#!/usr/bin/env python
##!/usr/bin/pvbatch --use-offscreen-rendering

################################################################################

try: paraview.simple
except: from paraview.simple import *
from ace3p_animation import Ace3pAnimataion, CameraPosition, CameraMovement
#except: 
#    import os
#    thisdir = __file__
#    thisdir = os.path.dirname(thisdir)
#    exec("from " + os.path.join(thisdir, "ace3p_animation") + 
#             " import Ace3pAnimataion, CameraPosition, CameraMovement")
from math import pi
################################################################################
# prevent goofy camera behavior
 
paraview.simple._DisableFirstRenderCameraReset()

################################################################################
anime = Ace3pAnimataion(mesh_file = '/home/jrr45/small_gap/small_gap_6.0.ncdf',
            mode_file = '/home/jrr45/small_gap/omega3p_results_6.0/omega3p.l0.m0000.1.4957339e+09.mod',
            trac3p_results = "/home/jrr45/small_gap/track3p_results_6.0/",
            reflections=['X'])

################################################################################
anime.load_files()

#find the center of the mesh
xbounds, ybounds, zbounds = anime.bounds
x = (xbounds[1]+xbounds[0])/2
y = (ybounds[1]+ybounds[0])/2
z = (zbounds[1]+zbounds[0])/2

position = CameraPosition(focus=[x, y, z], 
                          radius=1.0, theta=-pi/2, phi=pi/2, #theta=3*pi/2, phi=pi/2-pi/8, 
                          up=[-1.0,0.0,0.0], psi=0)

anime.play(camera_position=position, 
           camera_movements=[], GUI=False, antialias=2)
