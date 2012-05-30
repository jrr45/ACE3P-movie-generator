import sys, os, glob
from math import pi, sin, cos, sqrt
from paraview.simple import *
from numpy import cross, array

class CameraPosition(object):
    def __init__(self, focus, radius=0, theta=0, phi=0, up=[0.0,1.0,0.0], psi=0):
        self.view = GetRenderView()
        self.focus = focus
        self.radius = radius
        self.theta = theta
        self.phi = phi
        self.psi = psi
        self.tol = 10e-14
        
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
        print self.up
        
        
    def update_camera(self):
        print "updating"
        camera = GetActiveCamera()
        camera.SetFocalPoint(self.focus)
        camera.SetPosition(self.focus[0] + self.radius*sin(self.phi)*cos(self.theta), 
                           self.focus[1] + self.radius*sin(self.phi)*sin(self.theta),
                           self.focus[2] + self.radius*cos(self.phi))
        
        # rotation of the up vector around the camera vector
        [x, y, z] = self.up
        [u, v, w] = camera.GetPosition()
        calc1 = (u*x + v*y + w*z)*(cos(self.psi) - 1)
        calc2 = cos(self.psi)
        calc3 = sin(self.psi)
        camera.SetViewUp(-u*calc1 + x*calc2 + (-w*y + v*z)*calc3,
                         -v*calc1 + y*calc2 + ( w*x - u*z)*calc3,
                         -w*calc1 + z*calc2 + (-v*x + u*y)*calc3)
        print camera.GetViewUp()
        

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
        
        
class Ace3pAnimataion(object):
    """holds all information related to the animation of ACE3P"""
    def __init__(self, mesh_file,
             mode_file,
             trac3p_results="",
             reflections=[] # http://www.itk.org/Wiki/ParaView/Users_Guide/List_of_filters#Reflect
             ):
         
        # prevent goofy camera behavior
        paraview.simple._DisableFirstRenderCameraReset()
        
        # check the existance of the given files
        self.use_track3p = trac3p_results != ""
        self.mesh_file = mesh_file
        self.mode_file = mode_file
        self.trac3p_results = trac3p_results
        
        if not os.path.exists(mesh_file) or not os.path.isfile(mesh_file):
            print "mesh file was not found: %s" % mesh_file
            sys.exit(0)
        
        if not os.path.exists(mode_file) or not os.path.isfile(mode_file):
            print "mode file was not found: %s" % mode_file
            sys.exit(0)
            
        if self.use_track3p:
            particledir = os.path.join(trac3p_results, "PARTICLES")
            if not os.path.exists(particledir) or not os.path.isdir(particledir):
                print "trac3p results were not found (with any particles): %s" % particledir
                sys.exit(0)
        
            
            self.particleFilenames = glob.glob(os.path.join(particledir, "[0-9]*","partpath_ts*.ncdf"))
            self.particleFilenames.sort()
            
            if len(self.particleFilenames) == 0:
                print "No particles found at '%s'!" % particledir
                sys.exit(0)
                
        self.reflections = reflections
            
    def load_files(self):
        # Load mesh
        mesh = SLACDataReader(MeshFileName=self.mesh_file)
        mesh.ModeFileName = [self.mode_file]
        mesh.ReadInternalVolume = 1
        
        if self.use_track3p:
            # Load the particles
            self.p_ts = SLACParticleDataReader(FileName=self.particleFilenames)
            self.times = self.p_ts.TimestepValues
        else:
            self.times = mesh.TimestepValues.GetData()

        a3_efield_PVLookupTable = GetLookupTableForArray("efield", 3,
                                                          RGBPoints=[0.0, 0.0, 0.0, 1.0, 22.701986109850427, 1.0, 0.0, 0.0],
                                                          VectorMode='Magnitude',
                                                          ColorSpace='HSV',
                                                          ScalarRangeInitialized=1.0)
        a3_efield_PiecewiseFunction = CreatePiecewiseFunction()
        
        SetActiveSource(mesh)
        DataRepresentation1 = Show()
        DataRepresentation1.EdgeColor = [0.0, 0.0, 0.0]
        DataRepresentation1.ScalarOpacityUnitDistance = 0.034588492113178465
        DataRepresentation1.ScalarOpacityFunction = a3_efield_PiecewiseFunction
        DataRepresentation1.ColorArrayName = 'efield'
        DataRepresentation1.LookupTable = a3_efield_PVLookupTable
        DataRepresentation1.ColorAttributeType = 'POINT_DATA'
        DataRepresentation1.Representation = 'Wireframe'
        DataRepresentation1.BackfaceRepresentation = 'Surface'
        DataRepresentation1.Visibility = 0
        
        DataRepresentation2 = Show()
        DataRepresentation2.Visibility = 0
        DataRepresentation2.Representation = 'Outline'
        DataRepresentation2.EdgeColor = [0.0, 0.0, 0.0]
        
        #AnimationScene1 = GetAnimationScene()
        #AnimationScene1.EndTime = 6.2421777983841767e-07
        #AnimationScene1.PlayMode = 'Snap To TimeSteps'
        if self.use_track3p:
            SetActiveSource(self.p_ts)
            DataRepresentation3 = Show()
            DataRepresentation3.EdgeColor = [0.0, 0.0, 0.0]
            DataRepresentation3.Visibility = 0
            
        a3_efield_PVLookupTable.RGBPoints = [0.0, 0.0, 0.0, 1.0, 22.701986109850427, 1.0, 0.0, 0.0]
        
        
        SetActiveSource(mesh)
        ExtractBlock1 = ExtractBlock()
        ExtractBlock1.BlockIndices = [7]
        
        DataRepresentation4 = Show()
        DataRepresentation4.EdgeColor = [0.0, 0.0, 0.0]
        DataRepresentation4.ColorAttributeType = 'POINT_DATA'
        DataRepresentation4.ScalarOpacityFunction = a3_efield_PiecewiseFunction
        DataRepresentation4.ColorArrayName = 'efield'
        DataRepresentation4.ScalarOpacityUnitDistance = 0.037695299910287722
        DataRepresentation4.Texture = []
        DataRepresentation4.LookupTable = a3_efield_PVLookupTable
        DataRepresentation4.Visibility = 0
        
        # reflect mesh
        for plane in self.reflections:
            Reflect1 = Reflect()
            Reflect1.Plane = plane

        
        DataRepresentation5 = Show()
        DataRepresentation5.EdgeColor = [0.0, 0.0, 0.0]
        DataRepresentation5.ColorAttributeType = 'POINT_DATA'
        DataRepresentation5.ScalarOpacityFunction = a3_efield_PiecewiseFunction
        DataRepresentation5.ColorArrayName = 'efield'
        DataRepresentation5.ScalarOpacityUnitDistance = 0.030463782178688604
        DataRepresentation5.Texture = []
        DataRepresentation5.LookupTable = a3_efield_PVLookupTable
        DataRepresentation5.Representation = 'Wireframe'
        DataRepresentation5.BackfaceRepresentation = 'Surface'
        
        #reflect particles
        if self.use_track3p:
            SetActiveSource(self.p_ts)

            for plane in self.reflections:
                Reflect1 = Reflect()
                Reflect1.Plane = plane
            
            DataRepresentation7 = Show()
            #DataRepresentation7.ScalarOpacityUnitDistance = 0.043229431222641161
            DataRepresentation7.Texture = []
            DataRepresentation7.EdgeColor = [0.0, 0.0, 0.0]
    
        UpdatePipeline()
        bounds = mesh.GetDataInformation().GetBounds() 
        self.bounds = [[bounds[i*2], bounds[i*2+1]] for i in range(3)]
        
        # For some reason GetBounds does not include reflections so update bounds
        # with the new limits
        for plane in self.reflections:
            if plane == 'X':
                b = max(abs(self.bounds[0][0]), abs(self.bounds[0][1]))
                self.bounds[0] = [-b, b]
            elif plane == 'X min':
                self.bounds[0][0] -= self.bounds[0][1] - self.bounds[0][0]
            elif plane == 'X max':
                self.bounds[0][1] += self.bounds[0][1] - self.bounds[0][0]
            if plane == 'Y':
                b = max(abs(self.bounds[1][0]), abs(self.bounds[1][1]))
                self.bounds[1] = [-b, b]
            elif plane == 'Y min':
                self.bounds[1][0] -= self.bounds[1][1] - self.bounds[1][0]
            elif plane == 'Y max':
                self.bounds[1][1] += self.bounds[1][1] - self.bounds[1][0]
            if plane == 'Z':
                b = max(abs(self.bounds[2][0]), abs(self.bounds[2][1]))
                self.bounds[2] = [-b, b]
            elif plane == 'Z min':
                self.bounds[2][0] -= self.bounds[2][1] - self.bounds[2][0]
            elif plane == 'Z max':
                self.bounds[2][1] += self.bounds[2][1] - self.bounds[2][0]
        
    def play(self, camera_position, camera_movements=[], background_color=[0.0, 0.0, 0.0],
             antialias=1, peels=4, GUI=True, mono=True, stereo=False, speedup = 1):
        view = GetRenderView()
        view.Background = background_color # background color
        # Controls the number of passes in the rendering algorithm. A higher 
        # number of peels produces quality images, but increases rendering time.
        view.DepthPeeling = 0 
        view.MaximumNumberOfPeels = peels
        viewScale = 4
        view.ViewSize = [475, 300]
        
        if GUI:
            Render()
        else:
            times = self.times
            lentimes = len(times)
            digits = len(str(lentimes))
            
            if self.use_track3p:
                outFilenameEnd = "_particles.png"
            else:
                outFilenameEnd = "_mod.png"
        
            camera_position.update_camera()
            index = 0
            movelen = len(camera_movements)
               
            if movelen < index:
                movement = camera_movements[index]
            else:
                movement = None
            
            for ts in range(lentimes):
                if movement != None:
                    movement.update_camera(camera_position)
                    camera_position.update_camera()
                    movement.duration -= 1
                    
                    if movement.duration == 0:
                        if movelen < index:
                            movement = camera_movements[index]
                            index += 1
                        else:
                            movement = None

                if ts % speedup != 0:
                    continue
                
                
                now = times[ts] # simulation time
                outFilename = str(ts).zfill(digits) + outFilenameEnd
                
                print "===> Computing for '%s' ... of %i" % (outFilename, lentimes-1)
                
                ##########################################################################
                # update time to now, and show the results 
                view.ViewTime = now
        
                if stereo:
                    view.StereoRender = 1
        
                    view.StereoType = 'Left'
                    Show()
                    WriteImage("L_%s" % outFilename, Magnification=viewScale*antialias)
        
                    view.StereoType = 'Right'
                    Show()
                    WriteImage("R_%s" % outFilename, Magnification=viewScale*antialias)
        
                if mono:
                    view.StereoRender = 0
                    Show()
                    WriteImage(outFilename, Magnification=viewScale*antialias)
                

                
                
                # parametric fraction of animation seq.  in [0,1]
                #t = t_from_time(now) 
                
                ##########################################################################
                # move camera
                
                # position
                #(x, y, z) = xyz_from_ts(ts)
                #theta = (2.0 * pi / (48 * 30)) * ts
                #(x, z) = (x * cos(theta) + z * sin(theta), -x * sin(theta) + z * cos(theta))
                
                #x *= 1.4
                #y *= 1.4
                #z *= 1.4
        
                #view.CameraPosition = [x, y, z]
        
                # focal point
                #radius = 1.24
                #(fx, fy, fz) = getFocalPoint(x, y, z, radius)
                #(fx, fy, fz) = focus_from_ts(ts)
                #view.CameraFocalPoint = [fx, fy, fz]
        
                # compute up vector
                #(vx, vy, vz) = (fx - x, fy - y, fz - z) # view vector
                #(lx, ly, lz) = (0.0, 1.0, 0.0)    # "left" vector
        
              
                #ux = vy * lz - vz * ly # view x left = up
                #uy = vz * lx - vx * lz
                #uz = vx * ly - vy * lx
        
                #view.CameraViewUp = [ux, uy, uz]
        
                #print "ts:", ts
                #print "     x,  y,  z:", x, y, z
                #print "    fx, fy, fz:", fx, fy, fz
        
                #continue
