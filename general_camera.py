from math import sin, cos, sqrt
from numpy import cross, array
from camera import Camera
from paraview.simple import *

class CameraMovement(object):
    def __init__(self, duration, translation=[0,0,0], radius=0, theta=0, phi=0, psi=0):
        if duration < 1:
            raise ValueError("duration can not be less than 1")
        self.duration = duration
        self.translation = translation
        self.radius = radius
        self.theta = theta
        self.phi = phi
        self.psi = psi
        
    def update_camera(self, camera_position):
        camera_position.radius += self.radius/self.duration
        camera_position.theta  += self.theta /self.duration
        camera_position.phi    += self.phi   /self.duration
        camera_position.radius += self.radius/self.duration
        camera_position.radius += self.radius/self.duration
        
        #rotate the up vector - FIXME
        [u, v, w] = camera_position.up
        calc2 = w*(1-cos(self.theta /self.duration))
        calc3 = sin(self.theta /self.duration)
        camera_position.up = [u*calc2 + v*calc3,
                   v*calc2 - u*calc3,
                   w*calc2 + calc3]


class GeneralCamera(Camera):
    def __init__(self, *arguments, **keywords):
        """A general camera based on spherical coordinates with the following arguments:
           focus = [0.0, 0.0, 0.0]
           radius = 0
           theta = 0
           phi = 0
           up = [0.0, 1.0, 0.0]
           psi = 0
        """
        super( Camera, self ).__init__()

        self.view = GetRenderView()
        
        self.focus = [0, 0, 0]
        self.radius = 0
        self.theta = 0
        self.phi = 0
        self.psi = 0
        self.tol = 10e-14
        up = [0.0,1.0,0.0]
        self.camera_movements = []
        
        if "focus" in keywords:
            self.focus = keywords["focus"]
        if "camera_movement" in keywords:
            self.camera_movements = keywords["camera_movement"]
        if "radius" in keywords:
            self.radius = keywords["radius"]
        if "theta" in keywords:
            self.theta = keywords["theta"]
        if "phi" in keywords:
            self.phi = keywords["phi"]
        if "psi" in keywords:
            self.psi = keywords["psi"]
        if "up" in keywords:
            up = keywords["up"]
        
        
        #check to see if up vector is parallel to camera vector
        current = (self.focus[0] + self.radius*sin(self.phi)*cos(self.theta), 
                   self.focus[1] + self.radius*sin(self.phi)*sin(self.theta),
                   self.focus[2] + self.radius*cos(self.phi))
        vcamera = array([current])
        crossproduct = cross(array(up), vcamera)[0]
        if abs(crossproduct[0]) < self.tol and abs(crossproduct[1]) < self.tol and \
           abs(crossproduct[2]) < self.tol:
            print "up vector can not be parallel to initial camera focus"
            exit(-1)
        
        #normalize up vector
        norm = sqrt(up[0]**2 + up[1]**2 + up[2]**2)
        self.up = [up[0]/norm, up[1]/norm, up[2]/norm]
        
        
    def update_camera(self, framenumber):
        focus = self.focus
        theta = self.theta
        phi = self.phi
        radius = self.radius
        psi = 0
        
        # calculate sum of movement displacements
        for movement in self.camera_movements:
            # calculate the fraction of the movement to add to the position etc.
            delta_frames = 0
            if movement.duration > framenumber:
                delta_frames = movement.duration
            else:
                delta_frames = movement.duration
                
            if delta_frames == 0:
                break
            framenumber -= delta_frames
            
            fraction = delta_frames/movement.duration
        
            # update the total movement
            theta += movement.theta * fraction
            phi += movement.phi * fraction
            radius += movement.radius * fraction
            
            psi += movement.psi * fraction
            
            focus[0] += movement.translation[0] * fraction
            focus[1] += movement.translation[1] * fraction
            focus[2] += movement.translation[2] * fraction
        
        # finally set the camera at the new focus, position, and up vector 
        camera = GetActiveCamera()   
        camera.SetFocalPoint(focus)

        camera.SetPosition(focus[0] + radius*sin(phi)*cos(theta), 
                           focus[1] + radius*sin(phi)*sin(theta),
                           focus[2] + radius*cos(phi))
        
        # rotation of the up vector around the camera vector
        [x, y, z] = self.up
        [u, v, w] = camera.GetPosition()
        calc1 = (u*x + v*y + w*z)*(cos(self.psi) - 1)
        calc2 = cos(self.psi)
        calc3 = sin(self.psi)
        camera.SetViewUp(-u*calc1 + x*calc2 + (-w*y + v*z)*calc3,
                         -v*calc1 + y*calc2 + ( w*x - u*z)*calc3,
                         -w*calc1 + z*calc2 + (-v*x + u*y)*calc3)
        