# ParavView script to visualize the tumor. This script visualizes from t=0 to
# t=max and rotates once around the tumor during the visualization. The script
# must be used with the paraview python executable. The script in
# scripts/visualize-tumor-cells.sh wraps this script and can be used to
# visualize the tumor for multiple runs.
# Usage: pvpython paraview_tumor_slice.py <state_file> <transparent_background>
#                 <show_orientation_axes> <show_color_bar>
# transparent_background = 0  # 1 = transparent, 0 = grey
# show_orientation_axes = 1  # 1 = show, 0 = hide
# show_color_bar = True      # True = show, False = hide


# trace generated using paraview version 5.10.0
# import paraview
# paraview.compatibility.major = 5
# paraview.compatibility.minor = 10

#### import the simple module from the paraview
from paraview.simple import *
import os
import sys
import shutil


def visualize(
    filename, transparent_background, show_orientation_axes, show_color_bar
):
    # hard coded parameters
    output_folder = "slice"
    print(
        "<pvpython> Transparent background: {}".format(transparent_background)
    )
    print("<pvpython> Show orientation axes: {}".format(show_orientation_axes))
    print("<pvpython> Show color bar: {}".format(show_color_bar))

    # determine if we are running on an apple system
    is_apple = sys.platform == "darwin"

    #### disable automatic camera reset on 'Show'
    paraview.simple._DisableFirstRenderCameraReset()

    # get active view
    renderView1 = GetActiveViewOrCreate("RenderView")

    # destroy renderView1
    Delete(renderView1)
    del renderView1

    # Extract the folder name from filename
    folder = os.path.dirname(filename)

    # Extract the file name from filename
    file = os.path.basename(filename)

    filename = os.path.join(folder, file)
    print("<pvpython> Folder: " + folder)
    print("<pvpython> File: " + file)
    print("<pvpython> Loading state: " + filename)

    # load state
    LoadState(filename)
    print("<pvpython> Loading state done..")

    # find view
    renderView1 = FindViewOrCreate("RenderView1", viewtype="RenderView")

    # set active view
    SetActiveView(renderView1)

    # find source
    vessels = FindSource("Vessels")

    # hide data in view
    Hide(vessels, renderView1)

    # find source
    vEGFconcentration = FindSource("VEGF-concentration")

    # hide data in view
    Hide(vEGFconcentration, renderView1)

    # find source
    nutrientsconcentration = FindSource("Nutrients-concentration")

    # hide data in view
    Hide(nutrientsconcentration, renderView1)

    # get animation scene
    animationScene1 = GetAnimationScene()

    animationScene1.GoToLast()

    # find source
    tumorCells = FindSource("TumorCells")

    # set active source
    SetActiveSource(tumorCells)

    # get display properties
    tumorCellsDisplay = GetDisplayProperties(tumorCells, view=renderView1)
    renderView1.OrientationAxesVisibility = show_orientation_axes

    # set scalar coloring
    ColorBy(tumorCellsDisplay, ("POINTS", "cell_state_"))

    # rescale color and/or opacity maps used to include current data range
    tumorCellsDisplay.RescaleTransferFunctionToDataRange(True, False)

    # show color bar/color legend
    tumorCellsDisplay.SetScalarBarVisibility(renderView1, show_color_bar)

    # get color transfer function/color map for 'cell_state_'
    cell_state_LUT = GetColorTransferFunction("cell_state_")
    cell_state_LUT.AutomaticRescaleRangeMode = "Grow and update on 'Apply'"
    cell_state_LUT.InterpretValuesAsCategories = 0
    cell_state_LUT.AnnotationsInitialized = 0
    cell_state_LUT.ShowCategoricalColorsinDataRangeOnly = 0
    cell_state_LUT.RescaleOnVisibilityChange = 0
    cell_state_LUT.EnableOpacityMapping = 0
    cell_state_LUT.RGBPoints = [
        0.0,
        0.231373,
        0.298039,
        0.752941,
        2.0,
        0.865003,
        0.865003,
        0.865003,
        4.0,
        0.705882,
        0.0156863,
        0.14902,
    ]
    cell_state_LUT.UseLogScale = 0
    cell_state_LUT.UseOpacityControlPointsFreehandDrawing = 0
    cell_state_LUT.ShowDataHistogram = 0
    cell_state_LUT.AutomaticDataHistogramComputation = 0
    cell_state_LUT.DataHistogramNumberOfBins = 10
    cell_state_LUT.ColorSpace = "Diverging"
    cell_state_LUT.UseBelowRangeColor = 0
    cell_state_LUT.BelowRangeColor = [0.0, 0.0, 0.0]
    cell_state_LUT.UseAboveRangeColor = 0
    cell_state_LUT.AboveRangeColor = [0.5, 0.5, 0.5]
    cell_state_LUT.NanColor = [1.0, 1.0, 0.0]
    cell_state_LUT.NanOpacity = 1.0
    cell_state_LUT.Discretize = 1
    cell_state_LUT.NumberOfTableValues = 256
    cell_state_LUT.ScalarRangeInitialized = 1.0
    cell_state_LUT.HSVWrap = 0
    cell_state_LUT.VectorComponent = 0
    cell_state_LUT.VectorMode = "Magnitude"
    cell_state_LUT.AllowDuplicateScalars = 1
    cell_state_LUT.Annotations = []
    cell_state_LUT.ActiveAnnotatedValues = []
    cell_state_LUT.IndexedColors = []
    cell_state_LUT.IndexedOpacities = []

    # get opacity transfer function/opacity map for 'cell_state_'
    cell_state_PWF = GetOpacityTransferFunction("cell_state_")
    cell_state_PWF.Points = [0.0, 0.0, 0.5, 0.0, 4.0, 1.0, 0.5, 0.0]
    cell_state_PWF.AllowDuplicateScalars = 1
    cell_state_PWF.UseLogScale = 0
    cell_state_PWF.ScalarRangeInitialized = 1

    # Properties modified on cell_state_LUT
    cell_state_LUT.InterpretValuesAsCategories = 1
    cell_state_LUT.AnnotationsInitialized = 1

    # Properties modified on cell_state_LUT
    cell_state_LUT.Annotations = [
        "0",
        "Q",
        "1",
        "SG2",
        "2",
        "G1",
        "3",
        "H",
        "4",
        "D",
    ]

    # Properties modified on cell_state_LUT
    cell_state_LUT.IndexedColors = [
        0.803921568627451,
        0.5411764705882353,
        0.01568627450980392,
        0.050980392156862744,
        0.3843137254901961,
        0.0,
        0.6627450980392157,
        1.0,
        0.15294117647058825,
        0.396078431372549,
        0.396078431372549,
        0.40784313725490196,
        0.043137254901960784,
        0.043137254901960784,
        0.13725490196078433,
    ]
    cell_state_LUT.IndexedOpacities = [1.0, 1.0, 1.0, 1.0, 1.0]

    animationScene1.GoToFirst()

    # get animation track
    tumorCellsGlyphModeTrack = GetAnimationTrack(
        "GlyphMode", index=0, proxy=tumorCells
    )

    # create keyframes for this animation track

    # create a key frame
    keyFrame13404 = CompositeKeyFrame()
    keyFrame13404.KeyTime = 0.0
    keyFrame13404.KeyValues = [0.0]
    keyFrame13404.Interpolation = "Ramp"
    keyFrame13404.Base = 2.0
    keyFrame13404.StartPower = 0.0
    keyFrame13404.EndPower = 1.0
    keyFrame13404.Phase = 0.0
    keyFrame13404.Frequency = 1.0
    keyFrame13404.Offset = 0.0

    # create a key frame
    keyFrame13405 = CompositeKeyFrame()
    keyFrame13405.KeyTime = 1.0
    keyFrame13405.KeyValues = [0.0]
    keyFrame13405.Interpolation = "Ramp"
    keyFrame13405.Base = 2.0
    keyFrame13405.StartPower = 0.0
    keyFrame13405.EndPower = 1.0
    keyFrame13405.Phase = 0.0
    keyFrame13405.Frequency = 1.0
    keyFrame13405.Offset = 0.0

    # initialize the animation track
    tumorCellsGlyphModeTrack.TimeMode = "Normalized"
    tumorCellsGlyphModeTrack.StartTime = 0.0
    tumorCellsGlyphModeTrack.EndTime = 1.0
    tumorCellsGlyphModeTrack.Enabled = 1
    tumorCellsGlyphModeTrack.KeyFrames = [keyFrame13404, keyFrame13405]

    # get camera animation track for the view
    cameraAnimationCue1 = GetCameraTrack(view=renderView1)

    # create keyframes for this animation track

    # create a key frame
    keyFrame13410 = CameraKeyFrame()
    keyFrame13410.KeyTime = 0.0
    keyFrame13410.KeyValues = [0.0]
    keyFrame13410.Position = [
        -2.9318084716796875,
        0.21203231811523438,
        646.1783397197871,
    ]
    keyFrame13410.FocalPoint = [
        -2.9318084716796875,
        0.21203231811523438,
        1.5141716003417969,
    ]
    keyFrame13410.ViewUp = [0.0, 1.0, 0.0]
    keyFrame13410.ViewAngle = 30.0
    keyFrame13410.ParallelScale = 166.85136440448574
    # keyFrame13410.PositionPathPoints = [
    #     -2.9318084716796875,
    #     0.21203231811523438,
    #     502.7635803222656,
    #     291.69520169538924,
    #     0.21203231811523438,
    #     407.0334616767721,
    #     473.7847079823986,
    #     0.21203231811523438,
    #     156.40875731581025,
    #     473.7847079823985,
    #     0.21203231811523438,
    #     -153.38041411512648,
    #     291.6952016953892,
    #     0.21203231811523438,
    #     -404.0051184760882,
    #     -2.931808471679574,
    #     0.21203231811523438,
    #     -499.7352371215817,
    #     -297.5588186387483,
    #     0.21203231811523438,
    #     -404.00511847608834,
    #     -479.6483249257576,
    #     0.21203231811523438,
    #     -153.38041411512668,
    #     -479.64832492575766,
    #     0.21203231811523438,
    #     156.4087573158099,
    #     -297.55881863874845,
    #     0.21203231811523438,
    #     407.03346167677154,
    # ]
    keyFrame13410.FocalPathPoints = [
        -2.9318084716796875,
        0.21203231811523438,
        1.5141716003417969,
    ]
    keyFrame13410.PositionMode = "Path"
    keyFrame13410.FocalPointMode = "Path"
    keyFrame13410.ClosedFocalPath = 0
    keyFrame13410.ClosedPositionPath = 1

    keyFrame14203 = CameraKeyFrame()
    keyFrame14203.KeyTime = 0.0
    keyFrame14203.KeyValues = [0.0]
    keyFrame14203.Position = [-19.8793, -1.75063, 2837.85]
    keyFrame14203.FocalPoint = [-2.93181, -1.75063, -7.74627]
    keyFrame14203.ViewUp = [0.0, 1.0, 0.0]
    keyFrame14203.ViewAngle = 22.3773
    keyFrame14203.ParallelScale = 160.285
    keyFrame14203.PositionPathPoints = [
        -2.9318084716796875,
        0.21203231811523438,
        502.7635803222656,
        291.69520169538924,
        0.21203231811523438,
        407.0334616767721,
        473.7847079823986,
        0.21203231811523438,
        156.40875731581025,
        473.7847079823985,
        0.21203231811523438,
        -153.38041411512648,
        291.6952016953892,
        0.21203231811523438,
        -404.0051184760882,
        -2.931808471679574,
        0.21203231811523438,
        -499.7352371215817,
        -297.5588186387483,
        0.21203231811523438,
        -404.00511847608834,
        -479.6483249257576,
        0.21203231811523438,
        -153.38041411512668,
        -479.64832492575766,
        0.21203231811523438,
        156.4087573158099,
        -297.55881863874845,
        0.21203231811523438,
        407.03346167677154,
    ]
    keyFrame14203.FocalPathPoints = [
        -2.9318084716796875,
        0.21203231811523438,
        1.5141716003417969,
    ]
    keyFrame14203.PositionMode = "Path"
    keyFrame14203.FocalPointMode = "Path"
    keyFrame14203.ClosedFocalPath = 0
    keyFrame14203.ClosedPositionPath = 1

    # initialize the animation track
    cameraAnimationCue1.TimeMode = "Normalized"
    cameraAnimationCue1.StartTime = 0.0
    cameraAnimationCue1.EndTime = 1.0
    cameraAnimationCue1.Enabled = 1
    cameraAnimationCue1.Mode = "Interpolate Camera"
    cameraAnimationCue1.Interpolation = "Linear"
    cameraAnimationCue1.KeyFrames = [keyFrame13410, keyFrame14203]
    cameraAnimationCue1.DataSource = None

    animationScene1.GoToLast()

    # reset view to fit data bounds
    if is_apple:
        renderView1.ResetCamera(
            -385.5185852050781,
            470.3890686035156,
            -495.7220764160156,
            451.496826171875,
            -425.3491516113281,
            395.9783020019531,
            True,
        )
    else:
        renderView1.ResetCamera(
            -385.5185852050781,
            470.3890686035156,
            -495.7220764160156,
            451.496826171875,
            -425.3491516113281,
            395.9783020019531,
        )

    animationScene1.GoToPrevious()

    animationScene1.GoToFirst()

    # create a new 'Clip'
    clip1 = Clip(registrationName="Clip1", Input=tumorCells)
    clip1.ClipType = "Plane"
    clip1.HyperTreeGridClipper = "Plane"
    clip1.Scalars = ["POINTS", "cell_state_"]
    clip1.Value = 0.0
    clip1.Invert = 1
    clip1.Crinkleclip = 0
    clip1.Exact = 0

    # init the 'Plane' selected for 'ClipType'
    clip1.ClipType.Origin = [
        -2.9318084716796875,
        0.21203231811523438,
        1.5141716003417969,
    ]
    clip1.ClipType.Normal = [1.0, 0.0, 0.0]
    clip1.ClipType.Offset = 0.0

    # init the 'Plane' selected for 'HyperTreeGridClipper'
    clip1.HyperTreeGridClipper.Origin = [
        -2.9318084716796875,
        0.21203231811523438,
        1.5141716003417969,
    ]
    clip1.HyperTreeGridClipper.Normal = [1.0, 0.0, 0.0]
    clip1.HyperTreeGridClipper.Offset = 0.0

    # show data in view
    clip1Display = Show(clip1, renderView1, "UnstructuredGridRepresentation")

    # trace defaults for the display properties.
    clip1Display.Selection = None
    clip1Display.Representation = "Surface"
    clip1Display.ColorArrayName = ["POINTS", "cell_state_"]
    clip1Display.LookupTable = cell_state_LUT
    clip1Display.MapScalars = 1
    clip1Display.MultiComponentsMapping = 0
    clip1Display.InterpolateScalarsBeforeMapping = 1
    clip1Display.Opacity = 1.0
    clip1Display.PointSize = 2.0
    clip1Display.LineWidth = 1.0
    clip1Display.RenderLinesAsTubes = 0
    clip1Display.RenderPointsAsSpheres = 0
    clip1Display.Interpolation = "Gouraud"
    clip1Display.Specular = 0.0
    clip1Display.SpecularColor = [1.0, 1.0, 1.0]
    clip1Display.SpecularPower = 100.0
    clip1Display.Luminosity = 0.0
    clip1Display.Ambient = 0.0
    clip1Display.Diffuse = 1.0
    clip1Display.Roughness = 0.3
    clip1Display.Metallic = 0.0
    clip1Display.EdgeTint = [1.0, 1.0, 1.0]
    clip1Display.SelectTCoordArray = "None"
    clip1Display.SelectNormalArray = "Normals"
    clip1Display.SelectTangentArray = "None"
    clip1Display.Texture = None
    clip1Display.RepeatTextures = 1
    clip1Display.InterpolateTextures = 0
    clip1Display.SeamlessU = 0
    clip1Display.SeamlessV = 0
    clip1Display.UseMipmapTextures = 0
    clip1Display.BaseColorTexture = None
    clip1Display.NormalTexture = None
    clip1Display.NormalScale = 1.0
    clip1Display.MaterialTexture = None
    clip1Display.OcclusionStrength = 1.0
    clip1Display.EmissiveTexture = None
    clip1Display.EmissiveFactor = [1.0, 1.0, 1.0]
    clip1Display.FlipTextures = 0
    clip1Display.BackfaceRepresentation = "Follow Frontface"
    clip1Display.BackfaceAmbientColor = [1.0, 1.0, 1.0]
    clip1Display.BackfaceOpacity = 1.0
    clip1Display.Position = [0.0, 0.0, 0.0]
    clip1Display.Scale = [1.0, 1.0, 1.0]
    clip1Display.Orientation = [0.0, 0.0, 0.0]
    clip1Display.Origin = [0.0, 0.0, 0.0]
    clip1Display.CoordinateShiftScaleMethod = "Always Auto Shift Scale"
    clip1Display.Pickable = 1
    clip1Display.Triangulate = 0
    clip1Display.UseShaderReplacements = 0
    clip1Display.ShaderReplacements = ""
    clip1Display.NonlinearSubdivisionLevel = 1
    clip1Display.UseDataPartitions = 0
    clip1Display.OSPRayUseScaleArray = "All Approximate"
    clip1Display.OSPRayScaleArray = "Normals"
    clip1Display.OSPRayScaleFunction = "PiecewiseFunction"
    clip1Display.OSPRayMaterial = "None"
    clip1Display.Orient = 0
    clip1Display.OrientationMode = "Direction"
    clip1Display.SelectOrientationVectors = "None"
    clip1Display.Scaling = 0
    clip1Display.ScaleMode = "No Data Scaling Off"
    clip1Display.ScaleFactor = 20.049976348876953
    clip1Display.SelectScaleArray = "None"
    clip1Display.GlyphType = "Arrow"
    clip1Display.UseGlyphTable = 0
    clip1Display.GlyphTableIndexArray = "None"
    clip1Display.UseCompositeGlyphTable = 0
    clip1Display.UseGlyphCullingAndLOD = 0
    clip1Display.LODValues = []
    clip1Display.ColorByLODIndex = 0
    clip1Display.GaussianRadius = 1.0024988174438476
    clip1Display.ShaderPreset = "Sphere"
    clip1Display.CustomTriangleScale = 3
    clip1Display.CustomShader = """ // This custom shader code define a gaussian blur
    // Please take a look into vtkSMPointGaussianRepresentation.cxx
    // for other custom shader examples
    //VTK::Color::Impl
    float dist2 = dot(offsetVCVSOutput.xy,offsetVCVSOutput.xy);
    float gaussian = exp(-0.5*dist2);
    opacity = opacity*gaussian;
    """
    clip1Display.Emissive = 0
    clip1Display.ScaleByArray = 0
    clip1Display.SetScaleArray = ["POINTS", "Normals"]
    clip1Display.ScaleArrayComponent = "X"
    clip1Display.UseScaleFunction = 1
    clip1Display.ScaleTransferFunction = "PiecewiseFunction"
    clip1Display.OpacityByArray = 0
    clip1Display.OpacityArray = ["POINTS", "Normals"]
    clip1Display.OpacityArrayComponent = "X"
    clip1Display.OpacityTransferFunction = "PiecewiseFunction"
    clip1Display.DataAxesGrid = "GridAxesRepresentation"
    clip1Display.SelectionCellLabelBold = 0
    clip1Display.SelectionCellLabelColor = [0.0, 1.0, 0.0]
    clip1Display.SelectionCellLabelFontFamily = "Arial"
    clip1Display.SelectionCellLabelFontFile = ""
    clip1Display.SelectionCellLabelFontSize = 18
    clip1Display.SelectionCellLabelItalic = 0
    clip1Display.SelectionCellLabelJustification = "Left"
    clip1Display.SelectionCellLabelOpacity = 1.0
    clip1Display.SelectionCellLabelShadow = 0
    clip1Display.SelectionPointLabelBold = 0
    clip1Display.SelectionPointLabelColor = [1.0, 1.0, 0.0]
    clip1Display.SelectionPointLabelFontFamily = "Arial"
    clip1Display.SelectionPointLabelFontFile = ""
    clip1Display.SelectionPointLabelFontSize = 18
    clip1Display.SelectionPointLabelItalic = 0
    clip1Display.SelectionPointLabelJustification = "Left"
    clip1Display.SelectionPointLabelOpacity = 1.0
    clip1Display.SelectionPointLabelShadow = 0
    clip1Display.PolarAxes = "PolarAxesRepresentation"
    clip1Display.ScalarOpacityFunction = cell_state_PWF
    clip1Display.ScalarOpacityUnitDistance = 9.984539464260678
    clip1Display.UseSeparateOpacityArray = 0
    clip1Display.OpacityArrayName = ["POINTS", "Normals"]
    clip1Display.OpacityComponent = "X"
    clip1Display.SelectMapper = "Projected tetra"
    clip1Display.SamplingDimensions = [128, 128, 128]
    clip1Display.UseFloatingPointFrameBuffer = 1
    if is_apple:
        clip1Display.Anisotropy = 0.0
        clip1Display.AnisotropyRotation = 0.0
        clip1Display.BaseIOR = 1.5
        clip1Display.CoatStrength = 0.0
        clip1Display.CoatIOR = 2.0
        clip1Display.CoatRoughness = 0.0
        clip1Display.CoatColor = [1.0, 1.0, 1.0]
        clip1Display.ShowTexturesOnBackface = 1
        clip1Display.CoatNormalTexture = None
        clip1Display.CoatNormalScale = 1.0
        clip1Display.AnisotropyTexture = None
        clip1Display.BlockSelectors = ["/"]
        clip1Display.BlockColors = []
        clip1Display.BlockOpacities = []

    # init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
    clip1Display.OSPRayScaleFunction.Points = [
        0.0,
        0.0,
        0.5,
        0.0,
        1.0,
        1.0,
        0.5,
        0.0,
    ]
    clip1Display.OSPRayScaleFunction.UseLogScale = 0

    # init the 'Arrow' selected for 'GlyphType'
    clip1Display.GlyphType.TipResolution = 6
    clip1Display.GlyphType.TipRadius = 0.1
    clip1Display.GlyphType.TipLength = 0.35
    clip1Display.GlyphType.ShaftResolution = 6
    clip1Display.GlyphType.ShaftRadius = 0.03
    clip1Display.GlyphType.Invert = 0

    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    clip1Display.ScaleTransferFunction.Points = [
        -0.9749279618263245,
        0.0,
        0.5,
        0.0,
        0.9749279618263245,
        1.0,
        0.5,
        0.0,
    ]
    clip1Display.ScaleTransferFunction.UseLogScale = 0

    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    clip1Display.OpacityTransferFunction.Points = [
        -0.9749279618263245,
        0.0,
        0.5,
        0.0,
        0.9749279618263245,
        1.0,
        0.5,
        0.0,
    ]
    clip1Display.OpacityTransferFunction.UseLogScale = 0

    # init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
    clip1Display.DataAxesGrid.XTitle = "X Axis"
    clip1Display.DataAxesGrid.YTitle = "Y Axis"
    clip1Display.DataAxesGrid.ZTitle = "Z Axis"
    clip1Display.DataAxesGrid.XTitleFontFamily = "Arial"
    clip1Display.DataAxesGrid.XTitleFontFile = ""
    clip1Display.DataAxesGrid.XTitleBold = 0
    clip1Display.DataAxesGrid.XTitleItalic = 0
    clip1Display.DataAxesGrid.XTitleFontSize = 12
    clip1Display.DataAxesGrid.XTitleShadow = 0
    clip1Display.DataAxesGrid.XTitleOpacity = 1.0
    clip1Display.DataAxesGrid.YTitleFontFamily = "Arial"
    clip1Display.DataAxesGrid.YTitleFontFile = ""
    clip1Display.DataAxesGrid.YTitleBold = 0
    clip1Display.DataAxesGrid.YTitleItalic = 0
    clip1Display.DataAxesGrid.YTitleFontSize = 12
    clip1Display.DataAxesGrid.YTitleShadow = 0
    clip1Display.DataAxesGrid.YTitleOpacity = 1.0
    clip1Display.DataAxesGrid.ZTitleFontFamily = "Arial"
    clip1Display.DataAxesGrid.ZTitleFontFile = ""
    clip1Display.DataAxesGrid.ZTitleBold = 0
    clip1Display.DataAxesGrid.ZTitleItalic = 0
    clip1Display.DataAxesGrid.ZTitleFontSize = 12
    clip1Display.DataAxesGrid.ZTitleShadow = 0
    clip1Display.DataAxesGrid.ZTitleOpacity = 1.0
    clip1Display.DataAxesGrid.FacesToRender = 63
    clip1Display.DataAxesGrid.CullBackface = 0
    clip1Display.DataAxesGrid.CullFrontface = 1
    clip1Display.DataAxesGrid.ShowGrid = 0
    clip1Display.DataAxesGrid.ShowEdges = 1
    clip1Display.DataAxesGrid.ShowTicks = 1
    clip1Display.DataAxesGrid.LabelUniqueEdgesOnly = 1
    clip1Display.DataAxesGrid.AxesToLabel = 63
    clip1Display.DataAxesGrid.XLabelFontFamily = "Arial"
    clip1Display.DataAxesGrid.XLabelFontFile = ""
    clip1Display.DataAxesGrid.XLabelBold = 0
    clip1Display.DataAxesGrid.XLabelItalic = 0
    clip1Display.DataAxesGrid.XLabelFontSize = 12
    clip1Display.DataAxesGrid.XLabelShadow = 0
    clip1Display.DataAxesGrid.XLabelOpacity = 1.0
    clip1Display.DataAxesGrid.YLabelFontFamily = "Arial"
    clip1Display.DataAxesGrid.YLabelFontFile = ""
    clip1Display.DataAxesGrid.YLabelBold = 0
    clip1Display.DataAxesGrid.YLabelItalic = 0
    clip1Display.DataAxesGrid.YLabelFontSize = 12
    clip1Display.DataAxesGrid.YLabelShadow = 0
    clip1Display.DataAxesGrid.YLabelOpacity = 1.0
    clip1Display.DataAxesGrid.ZLabelFontFamily = "Arial"
    clip1Display.DataAxesGrid.ZLabelFontFile = ""
    clip1Display.DataAxesGrid.ZLabelBold = 0
    clip1Display.DataAxesGrid.ZLabelItalic = 0
    clip1Display.DataAxesGrid.ZLabelFontSize = 12
    clip1Display.DataAxesGrid.ZLabelShadow = 0
    clip1Display.DataAxesGrid.ZLabelOpacity = 1.0
    clip1Display.DataAxesGrid.XAxisNotation = "Mixed"
    clip1Display.DataAxesGrid.XAxisPrecision = 2
    clip1Display.DataAxesGrid.XAxisUseCustomLabels = 0
    clip1Display.DataAxesGrid.XAxisLabels = []
    clip1Display.DataAxesGrid.YAxisNotation = "Mixed"
    clip1Display.DataAxesGrid.YAxisPrecision = 2
    clip1Display.DataAxesGrid.YAxisUseCustomLabels = 0
    clip1Display.DataAxesGrid.YAxisLabels = []
    clip1Display.DataAxesGrid.ZAxisNotation = "Mixed"
    clip1Display.DataAxesGrid.ZAxisPrecision = 2
    clip1Display.DataAxesGrid.ZAxisUseCustomLabels = 0
    clip1Display.DataAxesGrid.ZAxisLabels = []
    clip1Display.DataAxesGrid.UseCustomBounds = 0
    clip1Display.DataAxesGrid.CustomBounds = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]

    # init the 'PolarAxesRepresentation' selected for 'PolarAxes'
    clip1Display.PolarAxes.Visibility = 0
    clip1Display.PolarAxes.Translation = [0.0, 0.0, 0.0]
    clip1Display.PolarAxes.Scale = [1.0, 1.0, 1.0]
    clip1Display.PolarAxes.Orientation = [0.0, 0.0, 0.0]
    clip1Display.PolarAxes.EnableCustomBounds = [0, 0, 0]
    clip1Display.PolarAxes.CustomBounds = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]
    clip1Display.PolarAxes.EnableCustomRange = 0
    clip1Display.PolarAxes.CustomRange = [0.0, 1.0]
    clip1Display.PolarAxes.PolarAxisVisibility = 1
    clip1Display.PolarAxes.RadialAxesVisibility = 1
    clip1Display.PolarAxes.DrawRadialGridlines = 1
    clip1Display.PolarAxes.PolarArcsVisibility = 1
    clip1Display.PolarAxes.DrawPolarArcsGridlines = 1
    clip1Display.PolarAxes.NumberOfRadialAxes = 0
    clip1Display.PolarAxes.AutoSubdividePolarAxis = 1
    clip1Display.PolarAxes.NumberOfPolarAxis = 0
    clip1Display.PolarAxes.MinimumRadius = 0.0
    clip1Display.PolarAxes.MinimumAngle = 0.0
    clip1Display.PolarAxes.MaximumAngle = 90.0
    clip1Display.PolarAxes.RadialAxesOriginToPolarAxis = 1
    clip1Display.PolarAxes.Ratio = 1.0
    clip1Display.PolarAxes.PolarAxisColor = [1.0, 1.0, 1.0]
    clip1Display.PolarAxes.PolarArcsColor = [1.0, 1.0, 1.0]
    clip1Display.PolarAxes.LastRadialAxisColor = [1.0, 1.0, 1.0]
    clip1Display.PolarAxes.SecondaryPolarArcsColor = [1.0, 1.0, 1.0]
    clip1Display.PolarAxes.SecondaryRadialAxesColor = [1.0, 1.0, 1.0]
    clip1Display.PolarAxes.PolarAxisTitleVisibility = 1
    clip1Display.PolarAxes.PolarAxisTitle = "Radial Distance"
    clip1Display.PolarAxes.PolarAxisTitleLocation = "Bottom"
    clip1Display.PolarAxes.PolarLabelVisibility = 1
    clip1Display.PolarAxes.PolarLabelFormat = "%-#6.3g"
    clip1Display.PolarAxes.PolarLabelExponentLocation = "Labels"
    clip1Display.PolarAxes.RadialLabelVisibility = 1
    clip1Display.PolarAxes.RadialLabelFormat = "%-#3.1f"
    clip1Display.PolarAxes.RadialLabelLocation = "Bottom"
    clip1Display.PolarAxes.RadialUnitsVisibility = 1
    clip1Display.PolarAxes.ScreenSize = 10.0
    clip1Display.PolarAxes.PolarAxisTitleOpacity = 1.0
    clip1Display.PolarAxes.PolarAxisTitleFontFamily = "Arial"
    clip1Display.PolarAxes.PolarAxisTitleFontFile = ""
    clip1Display.PolarAxes.PolarAxisTitleBold = 0
    clip1Display.PolarAxes.PolarAxisTitleItalic = 0
    clip1Display.PolarAxes.PolarAxisTitleShadow = 0
    clip1Display.PolarAxes.PolarAxisTitleFontSize = 12
    clip1Display.PolarAxes.PolarAxisLabelOpacity = 1.0
    clip1Display.PolarAxes.PolarAxisLabelFontFamily = "Arial"
    clip1Display.PolarAxes.PolarAxisLabelFontFile = ""
    clip1Display.PolarAxes.PolarAxisLabelBold = 0
    clip1Display.PolarAxes.PolarAxisLabelItalic = 0
    clip1Display.PolarAxes.PolarAxisLabelShadow = 0
    clip1Display.PolarAxes.PolarAxisLabelFontSize = 12
    clip1Display.PolarAxes.LastRadialAxisTextOpacity = 1.0
    clip1Display.PolarAxes.LastRadialAxisTextFontFamily = "Arial"
    clip1Display.PolarAxes.LastRadialAxisTextFontFile = ""
    clip1Display.PolarAxes.LastRadialAxisTextBold = 0
    clip1Display.PolarAxes.LastRadialAxisTextItalic = 0
    clip1Display.PolarAxes.LastRadialAxisTextShadow = 0
    clip1Display.PolarAxes.LastRadialAxisTextFontSize = 12
    clip1Display.PolarAxes.SecondaryRadialAxesTextOpacity = 1.0
    clip1Display.PolarAxes.SecondaryRadialAxesTextFontFamily = "Arial"
    clip1Display.PolarAxes.SecondaryRadialAxesTextFontFile = ""
    clip1Display.PolarAxes.SecondaryRadialAxesTextBold = 0
    clip1Display.PolarAxes.SecondaryRadialAxesTextItalic = 0
    clip1Display.PolarAxes.SecondaryRadialAxesTextShadow = 0
    clip1Display.PolarAxes.SecondaryRadialAxesTextFontSize = 12
    clip1Display.PolarAxes.EnableDistanceLOD = 1
    clip1Display.PolarAxes.DistanceLODThreshold = 0.7
    clip1Display.PolarAxes.EnableViewAngleLOD = 1
    clip1Display.PolarAxes.ViewAngleLODThreshold = 0.7
    clip1Display.PolarAxes.SmallestVisiblePolarAngle = 0.5
    clip1Display.PolarAxes.PolarTicksVisibility = 1
    clip1Display.PolarAxes.ArcTicksOriginToPolarAxis = 1
    clip1Display.PolarAxes.TickLocation = "Both"
    clip1Display.PolarAxes.AxisTickVisibility = 1
    clip1Display.PolarAxes.AxisMinorTickVisibility = 0
    clip1Display.PolarAxes.ArcTickVisibility = 1
    clip1Display.PolarAxes.ArcMinorTickVisibility = 0
    clip1Display.PolarAxes.DeltaAngleMajor = 10.0
    clip1Display.PolarAxes.DeltaAngleMinor = 5.0
    clip1Display.PolarAxes.PolarAxisMajorTickSize = 0.0
    clip1Display.PolarAxes.PolarAxisTickRatioSize = 0.3
    clip1Display.PolarAxes.PolarAxisMajorTickThickness = 1.0
    clip1Display.PolarAxes.PolarAxisTickRatioThickness = 0.5
    clip1Display.PolarAxes.LastRadialAxisMajorTickSize = 0.0
    clip1Display.PolarAxes.LastRadialAxisTickRatioSize = 0.3
    clip1Display.PolarAxes.LastRadialAxisMajorTickThickness = 1.0
    clip1Display.PolarAxes.LastRadialAxisTickRatioThickness = 0.5
    clip1Display.PolarAxes.ArcMajorTickSize = 0.0
    clip1Display.PolarAxes.ArcTickRatioSize = 0.3
    clip1Display.PolarAxes.ArcMajorTickThickness = 1.0
    clip1Display.PolarAxes.ArcTickRatioThickness = 0.5
    clip1Display.PolarAxes.Use2DMode = 0
    clip1Display.PolarAxes.UseLogAxis = 0

    # hide data in view
    Hide(tumorCells, renderView1)

    # show color bar/color legend
    clip1Display.SetScalarBarVisibility(renderView1, show_color_bar)

    # update the view to ensure updated data information
    renderView1.Update()

    # update the view to ensure updated data information
    renderView1.Update()

    # Properties modified on clip1.ClipType
    clip1.ClipType.Normal = [0.9999703645092188, 0.0, 0.0076987078980808085]

    # update the view to ensure updated data information
    renderView1.Update()

    # update the view to ensure updated data information
    renderView1.Update()

    # Properties modified on clip1.ClipType
    clip1.ClipType.Normal = [
        0.7103235705975832,
        0.2601608833213169,
        0.6540311459272964,
    ]

    # update the view to ensure updated data information
    renderView1.Update()

    # set active source
    SetActiveSource(tumorCells)

    # toggle 3D widget visibility (only when running from the GUI)
    Hide3DWidgets(proxy=clip1.ClipType)

    # update the view to ensure updated data information
    renderView1.Update()

    # create a new 'Clip'
    clip2 = Clip(registrationName="Clip2", Input=tumorCells)
    clip2.ClipType = "Plane"
    clip2.HyperTreeGridClipper = "Plane"
    clip2.Scalars = ["POINTS", "cell_state_"]
    clip2.Value = 0.0
    clip2.Invert = 1
    clip2.Crinkleclip = 0
    clip2.Exact = 0

    # init the 'Plane' selected for 'ClipType'
    clip2.ClipType.Origin = [
        -2.9318084716796875,
        0.21203231811523438,
        1.5141716003417969,
    ]
    clip2.ClipType.Normal = [1.0, 0.0, 0.0]
    clip2.ClipType.Offset = 0.0

    # init the 'Plane' selected for 'HyperTreeGridClipper'
    clip2.HyperTreeGridClipper.Origin = [
        -2.9318084716796875,
        0.21203231811523438,
        1.5141716003417969,
    ]
    clip2.HyperTreeGridClipper.Normal = [1.0, 0.0, 0.0]
    clip2.HyperTreeGridClipper.Offset = 0.0

    # show data in view
    clip2Display = Show(clip2, renderView1, "UnstructuredGridRepresentation")

    # trace defaults for the display properties.
    clip2Display.Selection = None
    clip2Display.Representation = "Surface"
    clip2Display.ColorArrayName = ["POINTS", "cell_state_"]
    clip2Display.LookupTable = cell_state_LUT
    clip2Display.MapScalars = 1
    clip2Display.MultiComponentsMapping = 0
    clip2Display.InterpolateScalarsBeforeMapping = 1
    clip2Display.Opacity = 1.0
    clip2Display.PointSize = 2.0
    clip2Display.LineWidth = 1.0
    clip2Display.RenderLinesAsTubes = 0
    clip2Display.RenderPointsAsSpheres = 0
    clip2Display.Interpolation = "Gouraud"
    clip2Display.Specular = 0.0
    clip2Display.SpecularColor = [1.0, 1.0, 1.0]
    clip2Display.SpecularPower = 100.0
    clip2Display.Luminosity = 0.0
    clip2Display.Ambient = 0.0
    clip2Display.Diffuse = 1.0
    clip2Display.Roughness = 0.3
    clip2Display.Metallic = 0.0
    clip2Display.EdgeTint = [1.0, 1.0, 1.0]
    clip2Display.SelectTCoordArray = "None"
    clip2Display.SelectNormalArray = "Normals"
    clip2Display.SelectTangentArray = "None"
    clip2Display.Texture = None
    clip2Display.RepeatTextures = 1
    clip2Display.InterpolateTextures = 0
    clip2Display.SeamlessU = 0
    clip2Display.SeamlessV = 0
    clip2Display.UseMipmapTextures = 0
    clip2Display.BaseColorTexture = None
    clip2Display.NormalTexture = None
    clip2Display.NormalScale = 1.0
    clip2Display.MaterialTexture = None
    clip2Display.OcclusionStrength = 1.0
    clip2Display.EmissiveTexture = None
    clip2Display.EmissiveFactor = [1.0, 1.0, 1.0]
    clip2Display.FlipTextures = 0
    clip2Display.BackfaceRepresentation = "Follow Frontface"
    clip2Display.BackfaceAmbientColor = [1.0, 1.0, 1.0]
    clip2Display.BackfaceOpacity = 1.0
    clip2Display.Position = [0.0, 0.0, 0.0]
    clip2Display.Scale = [1.0, 1.0, 1.0]
    clip2Display.Orientation = [0.0, 0.0, 0.0]
    clip2Display.Origin = [0.0, 0.0, 0.0]
    clip2Display.CoordinateShiftScaleMethod = "Always Auto Shift Scale"
    clip2Display.Pickable = 1
    clip2Display.Triangulate = 0
    clip2Display.UseShaderReplacements = 0
    clip2Display.ShaderReplacements = ""
    clip2Display.NonlinearSubdivisionLevel = 1
    clip2Display.UseDataPartitions = 0
    clip2Display.OSPRayUseScaleArray = "All Approximate"
    clip2Display.OSPRayScaleArray = "Normals"
    clip2Display.OSPRayScaleFunction = "PiecewiseFunction"
    clip2Display.OSPRayMaterial = "None"
    clip2Display.Orient = 0
    clip2Display.OrientationMode = "Direction"
    clip2Display.SelectOrientationVectors = "None"
    clip2Display.Scaling = 0
    clip2Display.ScaleMode = "No Data Scaling Off"
    clip2Display.ScaleFactor = 20.049976348876953
    clip2Display.SelectScaleArray = "None"
    clip2Display.GlyphType = "Arrow"
    clip2Display.UseGlyphTable = 0
    clip2Display.GlyphTableIndexArray = "None"
    clip2Display.UseCompositeGlyphTable = 0
    clip2Display.UseGlyphCullingAndLOD = 0
    clip2Display.LODValues = []
    clip2Display.ColorByLODIndex = 0
    clip2Display.GaussianRadius = 1.0024988174438476
    clip2Display.ShaderPreset = "Sphere"
    clip2Display.CustomTriangleScale = 3
    clip2Display.CustomShader = """ // This custom shader code define a gaussian blur
    // Please take a look into vtkSMPointGaussianRepresentation.cxx
    // for other custom shader examples
    //VTK::Color::Impl
    float dist2 = dot(offsetVCVSOutput.xy,offsetVCVSOutput.xy);
    float gaussian = exp(-0.5*dist2);
    opacity = opacity*gaussian;
    """
    clip2Display.Emissive = 0
    clip2Display.ScaleByArray = 0
    clip2Display.SetScaleArray = ["POINTS", "Normals"]
    clip2Display.ScaleArrayComponent = "X"
    clip2Display.UseScaleFunction = 1
    clip2Display.ScaleTransferFunction = "PiecewiseFunction"
    clip2Display.OpacityByArray = 0
    clip2Display.OpacityArray = ["POINTS", "Normals"]
    clip2Display.OpacityArrayComponent = "X"
    clip2Display.OpacityTransferFunction = "PiecewiseFunction"
    clip2Display.DataAxesGrid = "GridAxesRepresentation"
    clip2Display.SelectionCellLabelBold = 0
    clip2Display.SelectionCellLabelColor = [0.0, 1.0, 0.0]
    clip2Display.SelectionCellLabelFontFamily = "Arial"
    clip2Display.SelectionCellLabelFontFile = ""
    clip2Display.SelectionCellLabelFontSize = 18
    clip2Display.SelectionCellLabelItalic = 0
    clip2Display.SelectionCellLabelJustification = "Left"
    clip2Display.SelectionCellLabelOpacity = 1.0
    clip2Display.SelectionCellLabelShadow = 0
    clip2Display.SelectionPointLabelBold = 0
    clip2Display.SelectionPointLabelColor = [1.0, 1.0, 0.0]
    clip2Display.SelectionPointLabelFontFamily = "Arial"
    clip2Display.SelectionPointLabelFontFile = ""
    clip2Display.SelectionPointLabelFontSize = 18
    clip2Display.SelectionPointLabelItalic = 0
    clip2Display.SelectionPointLabelJustification = "Left"
    clip2Display.SelectionPointLabelOpacity = 1.0
    clip2Display.SelectionPointLabelShadow = 0
    clip2Display.PolarAxes = "PolarAxesRepresentation"
    clip2Display.ScalarOpacityFunction = cell_state_PWF
    clip2Display.ScalarOpacityUnitDistance = 9.984539464260678
    clip2Display.UseSeparateOpacityArray = 0
    clip2Display.OpacityArrayName = ["POINTS", "Normals"]
    clip2Display.OpacityComponent = "X"
    clip2Display.SelectMapper = "Projected tetra"
    clip2Display.SamplingDimensions = [128, 128, 128]
    clip2Display.UseFloatingPointFrameBuffer = 1
    if is_apple:
        clip2Display.Anisotropy = 0.0
        clip2Display.AnisotropyRotation = 0.0
        clip2Display.BaseIOR = 1.5
        clip2Display.CoatStrength = 0.0
        clip2Display.CoatIOR = 2.0
        clip2Display.CoatRoughness = 0.0
        clip2Display.CoatColor = [1.0, 1.0, 1.0]
        clip2Display.ShowTexturesOnBackface = 1
        clip2Display.CoatNormalTexture = None
        clip2Display.CoatNormalScale = 1.0
        clip2Display.AnisotropyTexture = None
        clip2Display.BlockSelectors = ["/"]
        clip2Display.BlockColors = []
        clip2Display.BlockOpacities = []

    # init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
    clip2Display.OSPRayScaleFunction.Points = [
        0.0,
        0.0,
        0.5,
        0.0,
        1.0,
        1.0,
        0.5,
        0.0,
    ]
    clip2Display.OSPRayScaleFunction.UseLogScale = 0

    # init the 'Arrow' selected for 'GlyphType'
    clip2Display.GlyphType.TipResolution = 6
    clip2Display.GlyphType.TipRadius = 0.1
    clip2Display.GlyphType.TipLength = 0.35
    clip2Display.GlyphType.ShaftResolution = 6
    clip2Display.GlyphType.ShaftRadius = 0.03
    clip2Display.GlyphType.Invert = 0

    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    clip2Display.ScaleTransferFunction.Points = [
        -0.9749279618263245,
        0.0,
        0.5,
        0.0,
        0.9749279618263245,
        1.0,
        0.5,
        0.0,
    ]
    clip2Display.ScaleTransferFunction.UseLogScale = 0

    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    clip2Display.OpacityTransferFunction.Points = [
        -0.9749279618263245,
        0.0,
        0.5,
        0.0,
        0.9749279618263245,
        1.0,
        0.5,
        0.0,
    ]
    clip2Display.OpacityTransferFunction.UseLogScale = 0

    # init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
    clip2Display.DataAxesGrid.XTitle = "X Axis"
    clip2Display.DataAxesGrid.YTitle = "Y Axis"
    clip2Display.DataAxesGrid.ZTitle = "Z Axis"
    clip2Display.DataAxesGrid.XTitleFontFamily = "Arial"
    clip2Display.DataAxesGrid.XTitleFontFile = ""
    clip2Display.DataAxesGrid.XTitleBold = 0
    clip2Display.DataAxesGrid.XTitleItalic = 0
    clip2Display.DataAxesGrid.XTitleFontSize = 12
    clip2Display.DataAxesGrid.XTitleShadow = 0
    clip2Display.DataAxesGrid.XTitleOpacity = 1.0
    clip2Display.DataAxesGrid.YTitleFontFamily = "Arial"
    clip2Display.DataAxesGrid.YTitleFontFile = ""
    clip2Display.DataAxesGrid.YTitleBold = 0
    clip2Display.DataAxesGrid.YTitleItalic = 0
    clip2Display.DataAxesGrid.YTitleFontSize = 12
    clip2Display.DataAxesGrid.YTitleShadow = 0
    clip2Display.DataAxesGrid.YTitleOpacity = 1.0
    clip2Display.DataAxesGrid.ZTitleFontFamily = "Arial"
    clip2Display.DataAxesGrid.ZTitleFontFile = ""
    clip2Display.DataAxesGrid.ZTitleBold = 0
    clip2Display.DataAxesGrid.ZTitleItalic = 0
    clip2Display.DataAxesGrid.ZTitleFontSize = 12
    clip2Display.DataAxesGrid.ZTitleShadow = 0
    clip2Display.DataAxesGrid.ZTitleOpacity = 1.0
    clip2Display.DataAxesGrid.FacesToRender = 63
    clip2Display.DataAxesGrid.CullBackface = 0
    clip2Display.DataAxesGrid.CullFrontface = 1
    clip2Display.DataAxesGrid.ShowGrid = 0
    clip2Display.DataAxesGrid.ShowEdges = 1
    clip2Display.DataAxesGrid.ShowTicks = 1
    clip2Display.DataAxesGrid.LabelUniqueEdgesOnly = 1
    clip2Display.DataAxesGrid.AxesToLabel = 63
    clip2Display.DataAxesGrid.XLabelFontFamily = "Arial"
    clip2Display.DataAxesGrid.XLabelFontFile = ""
    clip2Display.DataAxesGrid.XLabelBold = 0
    clip2Display.DataAxesGrid.XLabelItalic = 0
    clip2Display.DataAxesGrid.XLabelFontSize = 12
    clip2Display.DataAxesGrid.XLabelShadow = 0
    clip2Display.DataAxesGrid.XLabelOpacity = 1.0
    clip2Display.DataAxesGrid.YLabelFontFamily = "Arial"
    clip2Display.DataAxesGrid.YLabelFontFile = ""
    clip2Display.DataAxesGrid.YLabelBold = 0
    clip2Display.DataAxesGrid.YLabelItalic = 0
    clip2Display.DataAxesGrid.YLabelFontSize = 12
    clip2Display.DataAxesGrid.YLabelShadow = 0
    clip2Display.DataAxesGrid.YLabelOpacity = 1.0
    clip2Display.DataAxesGrid.ZLabelFontFamily = "Arial"
    clip2Display.DataAxesGrid.ZLabelFontFile = ""
    clip2Display.DataAxesGrid.ZLabelBold = 0
    clip2Display.DataAxesGrid.ZLabelItalic = 0
    clip2Display.DataAxesGrid.ZLabelFontSize = 12
    clip2Display.DataAxesGrid.ZLabelShadow = 0
    clip2Display.DataAxesGrid.ZLabelOpacity = 1.0
    clip2Display.DataAxesGrid.XAxisNotation = "Mixed"
    clip2Display.DataAxesGrid.XAxisPrecision = 2
    clip2Display.DataAxesGrid.XAxisUseCustomLabels = 0
    clip2Display.DataAxesGrid.XAxisLabels = []
    clip2Display.DataAxesGrid.YAxisNotation = "Mixed"
    clip2Display.DataAxesGrid.YAxisPrecision = 2
    clip2Display.DataAxesGrid.YAxisUseCustomLabels = 0
    clip2Display.DataAxesGrid.YAxisLabels = []
    clip2Display.DataAxesGrid.ZAxisNotation = "Mixed"
    clip2Display.DataAxesGrid.ZAxisPrecision = 2
    clip2Display.DataAxesGrid.ZAxisUseCustomLabels = 0
    clip2Display.DataAxesGrid.ZAxisLabels = []
    clip2Display.DataAxesGrid.UseCustomBounds = 0
    clip2Display.DataAxesGrid.CustomBounds = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]

    # init the 'PolarAxesRepresentation' selected for 'PolarAxes'
    clip2Display.PolarAxes.Visibility = 0
    clip2Display.PolarAxes.Translation = [0.0, 0.0, 0.0]
    clip2Display.PolarAxes.Scale = [1.0, 1.0, 1.0]
    clip2Display.PolarAxes.Orientation = [0.0, 0.0, 0.0]
    clip2Display.PolarAxes.EnableCustomBounds = [0, 0, 0]
    clip2Display.PolarAxes.CustomBounds = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]
    clip2Display.PolarAxes.EnableCustomRange = 0
    clip2Display.PolarAxes.CustomRange = [0.0, 1.0]
    clip2Display.PolarAxes.PolarAxisVisibility = 1
    clip2Display.PolarAxes.RadialAxesVisibility = 1
    clip2Display.PolarAxes.DrawRadialGridlines = 1
    clip2Display.PolarAxes.PolarArcsVisibility = 1
    clip2Display.PolarAxes.DrawPolarArcsGridlines = 1
    clip2Display.PolarAxes.NumberOfRadialAxes = 0
    clip2Display.PolarAxes.AutoSubdividePolarAxis = 1
    clip2Display.PolarAxes.NumberOfPolarAxis = 0
    clip2Display.PolarAxes.MinimumRadius = 0.0
    clip2Display.PolarAxes.MinimumAngle = 0.0
    clip2Display.PolarAxes.MaximumAngle = 90.0
    clip2Display.PolarAxes.RadialAxesOriginToPolarAxis = 1
    clip2Display.PolarAxes.Ratio = 1.0
    clip2Display.PolarAxes.PolarAxisColor = [1.0, 1.0, 1.0]
    clip2Display.PolarAxes.PolarArcsColor = [1.0, 1.0, 1.0]
    clip2Display.PolarAxes.LastRadialAxisColor = [1.0, 1.0, 1.0]
    clip2Display.PolarAxes.SecondaryPolarArcsColor = [1.0, 1.0, 1.0]
    clip2Display.PolarAxes.SecondaryRadialAxesColor = [1.0, 1.0, 1.0]
    clip2Display.PolarAxes.PolarAxisTitleVisibility = 1
    clip2Display.PolarAxes.PolarAxisTitle = "Radial Distance"
    clip2Display.PolarAxes.PolarAxisTitleLocation = "Bottom"
    clip2Display.PolarAxes.PolarLabelVisibility = 1
    clip2Display.PolarAxes.PolarLabelFormat = "%-#6.3g"
    clip2Display.PolarAxes.PolarLabelExponentLocation = "Labels"
    clip2Display.PolarAxes.RadialLabelVisibility = 1
    clip2Display.PolarAxes.RadialLabelFormat = "%-#3.1f"
    clip2Display.PolarAxes.RadialLabelLocation = "Bottom"
    clip2Display.PolarAxes.RadialUnitsVisibility = 1
    clip2Display.PolarAxes.ScreenSize = 10.0
    clip2Display.PolarAxes.PolarAxisTitleOpacity = 1.0
    clip2Display.PolarAxes.PolarAxisTitleFontFamily = "Arial"
    clip2Display.PolarAxes.PolarAxisTitleFontFile = ""
    clip2Display.PolarAxes.PolarAxisTitleBold = 0
    clip2Display.PolarAxes.PolarAxisTitleItalic = 0
    clip2Display.PolarAxes.PolarAxisTitleShadow = 0
    clip2Display.PolarAxes.PolarAxisTitleFontSize = 12
    clip2Display.PolarAxes.PolarAxisLabelOpacity = 1.0
    clip2Display.PolarAxes.PolarAxisLabelFontFamily = "Arial"
    clip2Display.PolarAxes.PolarAxisLabelFontFile = ""
    clip2Display.PolarAxes.PolarAxisLabelBold = 0
    clip2Display.PolarAxes.PolarAxisLabelItalic = 0
    clip2Display.PolarAxes.PolarAxisLabelShadow = 0
    clip2Display.PolarAxes.PolarAxisLabelFontSize = 12
    clip2Display.PolarAxes.LastRadialAxisTextOpacity = 1.0
    clip2Display.PolarAxes.LastRadialAxisTextFontFamily = "Arial"
    clip2Display.PolarAxes.LastRadialAxisTextFontFile = ""
    clip2Display.PolarAxes.LastRadialAxisTextBold = 0
    clip2Display.PolarAxes.LastRadialAxisTextItalic = 0
    clip2Display.PolarAxes.LastRadialAxisTextShadow = 0
    clip2Display.PolarAxes.LastRadialAxisTextFontSize = 12
    clip2Display.PolarAxes.SecondaryRadialAxesTextOpacity = 1.0
    clip2Display.PolarAxes.SecondaryRadialAxesTextFontFamily = "Arial"
    clip2Display.PolarAxes.SecondaryRadialAxesTextFontFile = ""
    clip2Display.PolarAxes.SecondaryRadialAxesTextBold = 0
    clip2Display.PolarAxes.SecondaryRadialAxesTextItalic = 0
    clip2Display.PolarAxes.SecondaryRadialAxesTextShadow = 0
    clip2Display.PolarAxes.SecondaryRadialAxesTextFontSize = 12
    clip2Display.PolarAxes.EnableDistanceLOD = 1
    clip2Display.PolarAxes.DistanceLODThreshold = 0.7
    clip2Display.PolarAxes.EnableViewAngleLOD = 1
    clip2Display.PolarAxes.ViewAngleLODThreshold = 0.7
    clip2Display.PolarAxes.SmallestVisiblePolarAngle = 0.5
    clip2Display.PolarAxes.PolarTicksVisibility = 1
    clip2Display.PolarAxes.ArcTicksOriginToPolarAxis = 1
    clip2Display.PolarAxes.TickLocation = "Both"
    clip2Display.PolarAxes.AxisTickVisibility = 1
    clip2Display.PolarAxes.AxisMinorTickVisibility = 0
    clip2Display.PolarAxes.ArcTickVisibility = 1
    clip2Display.PolarAxes.ArcMinorTickVisibility = 0
    clip2Display.PolarAxes.DeltaAngleMajor = 10.0
    clip2Display.PolarAxes.DeltaAngleMinor = 5.0
    clip2Display.PolarAxes.PolarAxisMajorTickSize = 0.0
    clip2Display.PolarAxes.PolarAxisTickRatioSize = 0.3
    clip2Display.PolarAxes.PolarAxisMajorTickThickness = 1.0
    clip2Display.PolarAxes.PolarAxisTickRatioThickness = 0.5
    clip2Display.PolarAxes.LastRadialAxisMajorTickSize = 0.0
    clip2Display.PolarAxes.LastRadialAxisTickRatioSize = 0.3
    clip2Display.PolarAxes.LastRadialAxisMajorTickThickness = 1.0
    clip2Display.PolarAxes.LastRadialAxisTickRatioThickness = 0.5
    clip2Display.PolarAxes.ArcMajorTickSize = 0.0
    clip2Display.PolarAxes.ArcTickRatioSize = 0.3
    clip2Display.PolarAxes.ArcMajorTickThickness = 1.0
    clip2Display.PolarAxes.ArcTickRatioThickness = 0.5
    clip2Display.PolarAxes.Use2DMode = 0
    clip2Display.PolarAxes.UseLogAxis = 0

    # hide data in view
    Hide(tumorCells, renderView1)

    # show color bar/color legend
    clip2Display.SetScalarBarVisibility(renderView1, show_color_bar)

    # update the view to ensure updated data information
    renderView1.Update()

    # update the view to ensure updated data information
    renderView1.Update()

    # Properties modified on clip2.ClipType
    clip2.ClipType.Normal = [
        -0.8524757964745938,
        0.08673548105853046,
        0.5155210691626022,
    ]

    # update the view to ensure updated data information
    renderView1.Update()

    # set active source
    SetActiveSource(tumorCells)

    # toggle 3D widget visibility (only when running from the GUI)
    Hide3DWidgets(proxy=clip2.ClipType)

    # update the view to ensure updated data information
    renderView1.Update()

    # get layout
    layout1 = GetLayout()

    # layout/tab size in pixels
    layout1.SetSize(1333, 942)

    # current camera placement for renderView1
    renderView1.CameraPosition = [-2.93181, 0.212032, 646.178]
    renderView1.CameraFocalPoint = [-2.93181, 0.212032, 1.51417]
    renderView1.CameraParallelScale = 166.851

    # Enable OSPRay for rendering on server
    if not is_apple:
        pm = paraview.servermanager.vtkSMProxyManager
        if pm.GetVersionMajor() == 5 and pm.GetVersionMinor() < 7:
            renderView1.EnableOSPRay = 1
            renderView1.OSPRayRenderer = 'pathtracer'
        else:
            renderView1.EnableRayTracing = 1
            renderView1.BackEnd = 'OSPRay pathtracer'
            renderView1.Denoise = 1
        # Properties modified on renderView1
        renderView1.Shadows = 1
        # Properties modified on renderView1
        renderView1.SamplesPerPixel = 20

    # save animation
    print("<pvpython> Create folder ..")
    output_folder = output_folder + "_bg{}_cb{}_ax{}".format(
        transparent_background, int(show_color_bar), show_orientation_axes
    )
    animation_folder = os.path.join(folder, output_folder)
    if not os.path.exists(animation_folder):
        os.makedirs(animation_folder)
    else:
        print("<pvpython> Folder already exists ..")
        print("<pvpython> Delete folder ..")
        shutil.rmtree(animation_folder)
        print("<pvpython> Create folder ..")
        os.makedirs(animation_folder)
    print("<pvpython> Saving animation ..")
    animation_path = os.path.join(folder, output_folder, "img.png")
    SaveAnimation(
        animation_path,
        renderView1,
        ImageResolution=[2704, 1520],
        FontScaling="Scale fonts proportionally",
        OverrideColorPalette="",
        StereoMode="No change",
        TransparentBackground=transparent_background,
        FrameRate=1,
        FrameWindow=[0, 99],
        # PNG options
        CompressionLevel="1",
        SuffixFormat=".%04d",
    )


def main(argc, argv):
    if argc != 5:
        print(
            "Usage: visualize.py <state_file> <transparent_background> "
            + "<show_orientation_axes> <show_color_bar>"
        )
        return
    filename = argv[1]
    transparent_background = int(argv[2])
    show_orientation_axes = int(argv[3])
    show_color_bar = bool(int(argv[4]))
    # check if file filename exists
    if not os.path.isfile(filename):
        print("File {} does not exist".format(filename))
        return
    visualize(
        filename, transparent_background, show_orientation_axes, show_color_bar
    )


if __name__ == "__main__":
    main(len(sys.argv), sys.argv)
