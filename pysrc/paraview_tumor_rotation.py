# ParavView script to visualize the tumor. This script visualizes from t=0 to
# t=max and rotates once around the tumor during the visualization. The script
# must be used with the paraview python executable. The script in
# scripts/visualize-tumor-cells.sh wraps this script and can be used to
# visualize the tumor for multiple runs.
# Usage: pvpython paraview_tumor_rotation.py <state_file> <transparent_background>
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
    output_folder = "rotation"
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

    # reset view to fit data
    if is_apple:
        renderView1.ResetCamera(True)

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
    cell_state_LUT.IndexedOpacities = [1.0, 1.0, 1.0, 1.0, 1.0]

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
        42.43524169921875,
        -22.112625122070312,
        2917.9012382262445,
    ]
    keyFrame13410.FocalPoint = [
        42.43524169921875,
        -22.112625122070312,
        -14.6854248046875,
    ]
    keyFrame13410.ViewUp = [0.0, 1.0, 0.0]
    keyFrame13410.ViewAngle = 23.354564755838638
    keyFrame13410.ParallelScale = 759.0092798060535
    keyFrame13410.PositionPathPoints = [
        42.4352,
        -22.1126,
        2917.9,
        2335.2227907461,
        -22.1126,
        1813.751689979815,
        2901.4945613168984,
        -22.1126,
        -667.2470421146517,
        1314.8363186335603,
        -22.1126,
        -2656.8535478651634,
        -1229.9659186335598,
        -22.1126,
        -2656.853547865164,
        -2816.624161316899,
        -22.1126,
        -667.2470421146527,
        -2250.3523907461017,
        -22.1126,
        1813.751689979815,
    ]
    keyFrame13410.FocalPathPoints = [42.4352, -22.1126, -14.6854]
    keyFrame13410.PositionMode = "Path"
    keyFrame13410.FocalPointMode = "Path"
    keyFrame13410.ClosedFocalPath = 0
    keyFrame13410.ClosedPositionPath = 1

    # create a key frame
    keyFrame13411 = CameraKeyFrame()
    keyFrame13411.KeyTime = 1.0
    keyFrame13411.KeyValues = [0.0]
    keyFrame13411.Position = [
        42.43524169921875,
        -22.112625122070312,
        2917.9012382262445,
    ]
    keyFrame13411.FocalPoint = [
        42.43524169921875,
        -22.112625122070312,
        -14.6854248046875,
    ]
    keyFrame13411.ViewUp = [0.0, 1.0, 0.0]
    keyFrame13411.ViewAngle = 23.354564755838638
    keyFrame13411.ParallelScale = 759.0092798060535
    keyFrame13411.PositionPathPoints = [
        5.0,
        0.0,
        0.0,
        5.0,
        5.0,
        0.0,
        5.0,
        0.0,
        0.0,
    ]
    keyFrame13411.FocalPathPoints = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0]
    keyFrame13411.PositionMode = "Path"
    keyFrame13411.FocalPointMode = "Path"
    keyFrame13411.ClosedFocalPath = 0
    keyFrame13411.ClosedPositionPath = 0

    # initialize the animation track
    cameraAnimationCue1.TimeMode = "Normalized"
    cameraAnimationCue1.StartTime = 0.0
    cameraAnimationCue1.EndTime = 1.0
    cameraAnimationCue1.Enabled = 1
    cameraAnimationCue1.Mode = "Path-based"
    cameraAnimationCue1.Interpolation = "Spline"
    cameraAnimationCue1.KeyFrames = [keyFrame13410, keyFrame13411]
    cameraAnimationCue1.DataSource = None

    # get layout
    layout1 = GetLayout()

    # layout/tab size in pixels
    layout1.SetSize(1333, 942)

    # current camera placement for renderView1
    renderView1.CameraPosition = [
        42.43524169921875,
        -22.112625122070312,
        2917.9012382262445,
    ]
    renderView1.CameraFocalPoint = [
        42.43524169921875,
        -22.112625122070312,
        -14.6854248046875,
    ]
    renderView1.CameraViewAngle = 23.354564755838638
    renderView1.CameraParallelScale = 759.0092798060535

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

    # # ================================================================
    # # addendum: following script captures some of the application
    # # state to faithfully reproduce the visualization during playback
    # # ================================================================

    # # --------------------------------
    # # saving layout sizes for layouts

    # # layout/tab size in pixels
    # layout1.SetSize(935, 822)

    # # -----------------------------------
    # # saving camera placements for views

    # # current camera placement for renderView1
    # renderView1.CameraPosition = [42.4352, -22.1126, 2917.9]
    # renderView1.CameraFocalPoint = [42.4352, -22.1126, -14.6854]
    # renderView1.CameraViewAngle = 23.450854700854702
    # renderView1.CameraParallelScale = 759.0092798060535

    # # --------------------------------------------
    # # uncomment the following to render all views
    # # RenderAllViews()
    # # alternatively, if you want to write images, you can use SaveScreenshot(...).


def main(argc, argv):
    if argc != 5:
        print(
            "Usage: visualize.py <filename> <transparent_background> "
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
