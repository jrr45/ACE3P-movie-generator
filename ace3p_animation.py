import sys, os, glob
from paraview.simple import *
        
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

        # Load mesh
        self.mesh = SLACDataReader(MeshFileName=self.mesh_file)
        self.mesh.ModeFileName = [self.mode_file]
        self.mesh.ReadInternalVolume = 1
        
        if self.use_track3p:
            # Load the particles
            self.p_ts = SLACParticleDataReader(FileName=self.particleFilenames)
            self.times = self.p_ts.TimestepValues
        else:
            self.times = self.mesh.TimestepValues.GetData()

        a3_efield_PVLookupTable = GetLookupTableForArray("efield", 3,
                                                          RGBPoints=[0.0, 0.0, 0.0, 1.0, 22.701986109850427, 1.0, 0.0, 0.0],
                                                          VectorMode='Magnitude',
                                                          ColorSpace='HSV',
                                                          ScalarRangeInitialized=1.0)
        a3_efield_PiecewiseFunction = CreatePiecewiseFunction()
        
        SetActiveSource(self.mesh)
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
        
        #outline box
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
        
        SetActiveSource(self.mesh)
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
            DataRepresentation7.ScalarOpacityUnitDistance = 0.043229431222641161
            DataRepresentation7.Texture = []
            DataRepresentation7.EdgeColor = [0.0, 0.0, 0.0]
    
        UpdatePipeline()
        
    def get_bounds(self):
        bounds = self.mesh.GetDataInformation().GetBounds() 
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
        return self.bounds
                
    def get_center(self):
        xbounds, ybounds, zbounds = self.get_bounds()
        x = (xbounds[1]+xbounds[0])/2
        y = (ybounds[1]+ybounds[0])/2
        z = (zbounds[1]+zbounds[0])/2
        return x, y, z
        
    def play(self, camera, background_color=[0.0, 0.0, 0.0],
             antialias=1, peels=4, GUI=True, mono=True, stereo=False, speedup = 1):
        view = GetRenderView()
        view.Background = background_color # background color
        # Controls the number of passes in the rendering algorithm. A higher 
        # number of peels produces quality images, but increases rendering time.
        view.DepthPeeling = 0 
        view.MaximumNumberOfPeels = peels
        viewScale = 4
        view.ViewSize = [480, 270]
        
        if GUI:
            camera.update_camera(0)
            Render()
        else:
            times = self.times
            lentimes = len(times)
            digits = len(str(lentimes))
            
            if self.use_track3p:
                outFilenameEnd = "_particles.png"
            else:
                outFilenameEnd = "_mod.png"
                
            if not os.path.exists("movie_images"):
                os.mkdir("movie_images")
            
            for ts in range(lentimes):

                if ts % speedup != 0:
                    continue
                
                camera.update_camera(ts)
                
                now = times[ts] # simulation time
                outFilename = "movie_images/" + str(ts).zfill(digits) + outFilenameEnd
                
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
                    
            if stereo:
                print "boo"
            else:
                print "now run:"
                print ("ffmpeg -i movie_images/%0" + str(digits) + "d_particles.png -s 1920x1080 " +
                      "movie.mp4")

