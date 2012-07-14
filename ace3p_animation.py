import sys, os, glob
from paraview.simple import *
        
class Ace3pAnimataion(object):
    """holds all information related to the animation of ACE3P"""
    def __init__(self, mesh_file,
             mode_files=[],
             trac3p_results="",
             reflections=[], # http://www.itk.org/Wiki/ParaView/Users_Guide/List_of_filters#Reflect
             sidesetIDs=[6],
             field='efield',
             VectorMode='Magnitude',#['Magnitude', 'Component']
             VectorComponent='X',#['X','Y', 'Z'] only used if VectorMode='Component'
             ColorSpace='HSV',
             RGBmin=[0.0, 0.0, 1.0],
             RGBmax=[1.0, 0.0, 0.0],
             Representation='Wireframe', #['Outline', 'Points', 'Wireframe', 'Surface', 'Surface With Edges'] 
             BackfaceRepresentation='Surface'
             ):
         
        # prevent goofy camera behavior
        paraview.simple._DisableFirstRenderCameraReset()
        
        # check the existance of the given files
        mode_files.sort()
        self.use_track3p = trac3p_results != ""
        self.mesh_file = mesh_file
        self.mode_files = mode_files
        self.trac3p_results = trac3p_results
        
        if not os.path.exists(mesh_file) or not os.path.isfile(mesh_file):
            print "mesh file was not found: %s" % mesh_file
            sys.exit(0)
        
        if len(mode_files) == 0:
            print "there must be at least one mode file"
            sys.exit(0)
                
        for mode_file in mode_files:
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
        self.mesh.ModeFileName = self.mode_files
        self.mesh.ReadInternalVolume = 1
        self.mesh.UpdatePipelineInformation()
        
        if self.use_track3p:
            # Load the particles
            self.p_ts = SLACParticleDataReader(FileName=self.particleFilenames)
            
        field_data = self.mesh.PointData.GetArray(field)
        if field_data is None:
            print "invalid field name: %s" % field
            sys.exit(0)
        
        # set of the coloring of the mesh
        field_range = field_data.GetRange()
        RGBPoints=[max(0, field_range[0]), # first point is the min of the range
                   RGBmin[0], RGBmin[1], RGBmin[2], #the color at the minimum
                   max(abs(field_range[0]), abs(field_range[1])), #max of the field range
                   RGBmax[0], RGBmax[1], RGBmax[2]] #the color at the maximum

        if VectorMode == 'Component':
            comp = {'X': 0, 'Y': 1, 'Z': 2, 'x': 0, 'y': 1, 'z': 2}
            field_PVLookupTable = GetLookupTableForArray(field, 3,
                                                         RGBPoints=RGBPoints,
                                                         VectorMode='Component',
                                                         VectorComponent=comp[VectorComponent],
                                                         ColorSpace=ColorSpace,
                                                         ScalarRangeInitialized=1.0)
        else:
            field_PVLookupTable = GetLookupTableForArray(field, 3,
                                                         RGBPoints=RGBPoints,
                                                         VectorMode=VectorMode,
                                                         ColorSpace=ColorSpace,
                                                         ScalarRangeInitialized=1.0)
        field_PiecewiseFunction = CreatePiecewiseFunction()
        
        # select the mesh surfaces that match the given Cubit sidesets
        SetActiveSource(self.mesh)
        ExtractBlock1 = ExtractBlock()
        ExtractBlock1.BlockIndices = [x+1 for x in sidesetIDs]
        
        # reflect mesh oven the given planes
        for plane in self.reflections:
            Reflect1 = Reflect()
            Reflect1.Plane = plane
            
        # magic number???
        ScalarOpacityUnitDistance = 1.0
        
        # Draw the fields
        DataRepresentation5 = Show()
        DataRepresentation5.EdgeColor = [0.0, 0.0, 0.0]
        DataRepresentation5.ColorAttributeType = 'POINT_DATA'
        DataRepresentation5.ScalarOpacityFunction = field_PiecewiseFunction
        DataRepresentation5.ColorArrayName = field
        DataRepresentation5.ScalarOpacityUnitDistance = ScalarOpacityUnitDistance
        DataRepresentation5.Texture = []
        DataRepresentation5.LookupTable = field_PVLookupTable
        DataRepresentation5.Representation = Representation
        DataRepresentation5.BackfaceRepresentation = BackfaceRepresentation
        
        #reflect particles
        if self.use_track3p:
            SetActiveSource(self.p_ts)

            for plane in self.reflections:
                Reflect1 = Reflect()
                Reflect1.Plane = plane
            
            DataRepresentation7 = Show()
            DataRepresentation7.ScalarOpacityUnitDistance = ScalarOpacityUnitDistance
            DataRepresentation7.Texture = []
            DataRepresentation7.EdgeColor = [0.0, 0.0, 0.0]
    
        UpdatePipeline()
        
        #get time information
        if self.use_track3p:
            self.times = self.p_ts.TimestepValues
        else:
            SetActiveSource(self.mesh)
            self.times = self.mesh.TimestepValues
            # Some mode files will not have TimeSteps so generate some
            if len(self.times) == 0:
                timerange = self.mesh.GetProperty('TimeRange') 
                steps = 100
                self.times = [timerange[1]*i/steps + timerange[0] for i in range(steps)]
                
        
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
    
    def get_frame_range(self):
        return [0, len(self.times)]
        
    def play(self, camera, background_color=[0.0, 0.0, 0.0],
             antialias=1, peels=4, GUI=True, mono=True, stereo=False, 
             speedup = 1, frame_range=None):
        view = GetRenderView()
        view.Background = background_color # background color
        # Controls the number of passes in the rendering algorithm. A higher 
        # number of peels produces quality images, but increases rendering time.
        view.DepthPeeling = 0 
        view.MaximumNumberOfPeels = peels
        viewScale = 4
        view.ViewSize = [480, 270]
        view.UseLight = 1
        
        if GUI:
            camera.update_camera(0)
            Render()
            return
        
        times = self.times
        lentimes = len(times)
        digits = len(str(lentimes))
        
        if self.use_track3p:
            outFilenameEnd = "_particles.png"
        else:
            outFilenameEnd = "_mod.png"
            
        if not os.path.exists("movie_images"):
            os.mkdir("movie_images")
        
        # set the range to be used in the loop
        if frame_range is not None:
            if frame_range[0] >= 0 and frame_range[1] <= lentimes and \
               frame_range[0] < frame_range[1]:
                timerange = range(frame_range[0], frame_range[1])
            else:
                raise ValueError("invalid frame range")
        else:
            timerange = range(lentimes)
        
        for ts in timerange:

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
            print ("The ffmpeg command will vary depending on how you are constructing the " +
                  "3D movie")
            print "the command for a 2D movie would be:"
        else:
            print "now run:"
        print ("ffmpeg -i movie_images/%0" + str(digits) + outFilenameEnd + " -s 1920x1080 " +
               "movie.mp4")