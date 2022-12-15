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
    filename,
    transparent_background,
    show_orientation_axes,
    show_color_bar,
    overwrite,
):
    print(
        "<pvpython> Transparent background: {}".format(transparent_background)
    )
    print("<pvpython> Show orientation axes: {}".format(show_orientation_axes))
    print("<pvpython> Show color bar: {}".format(show_color_bar))
    print("<pvpython> Overwrite: {}".format(overwrite))

    # determine if we are running on an apple system
    is_apple = sys.platform == "darwin"

    #### disable automatic camera reset on 'Show'
    paraview.simple._DisableFirstRenderCameraReset()

    # get active view
    renderView1 = GetActiveViewOrCreate("RenderView")

    # destroy renderView1
    Delete(renderView1)
    del renderView1

    # load state
    LoadState(filename)

    # find view
    renderView1 = FindViewOrCreate("RenderView1", viewtype="RenderView")

    # set active view
    SetActiveView(renderView1)

    # find source
    tumorCells = FindSource("TumorCells")

    # hide data in view
    Hide(tumorCells, renderView1)

    # find source
    nutrientsconcentration = FindSource("Nutrients-concentration")

    # hide data in view
    Hide(nutrientsconcentration, renderView1)

    # find source
    dOXconcentration = FindSource("DOX-concentration")

    # hide data in view
    Hide(dOXconcentration, renderView1)

    # find source
    tRAconcentration = FindSource("TRA-concentration")

    # hide data in view
    Hide(tRAconcentration, renderView1)

    # find source
    vEGFconcentration = FindSource("VEGF-concentration")

    # hide data in view
    Hide(vEGFconcentration, renderView1)

    # set active source
    SetActiveSource(vEGFconcentration)

    # get display properties
    vEGFconcentrationDisplay = GetDisplayProperties(
        vEGFconcentration, view=renderView1
    )

    # toggle 3D widget visibility (only when running from the GUI)
    Show3DWidgets(proxy=vEGFconcentrationDisplay)

    # toggle 3D widget visibility (only when running from the GUI)
    Hide3DWidgets(proxy=vEGFconcentrationDisplay.SliceFunction)

    # toggle 3D widget visibility (only when running from the GUI)
    Hide3DWidgets(proxy=vEGFconcentrationDisplay)

    # show data in view
    vEGFconcentrationDisplay = Show(
        vEGFconcentration, renderView1, "UniformGridRepresentation"
    )

    # show color bar/color legend
    vEGFconcentrationDisplay.SetScalarBarVisibility(renderView1, show_color_bar)

    # get color transfer function/color map for 'SubstanceConcentration'
    substanceConcentrationLUT = GetColorTransferFunction(
        "SubstanceConcentration"
    )
    substanceConcentrationLUT.AutomaticRescaleRangeMode = (
        "Grow and update on 'Apply'"
    )
    substanceConcentrationLUT.InterpretValuesAsCategories = 0
    substanceConcentrationLUT.AnnotationsInitialized = 0
    substanceConcentrationLUT.ShowCategoricalColorsinDataRangeOnly = 0
    substanceConcentrationLUT.RescaleOnVisibilityChange = 0
    substanceConcentrationLUT.EnableOpacityMapping = 0
    substanceConcentrationLUT.RGBPoints = [
        0.0,
        0.231373,
        0.298039,
        0.752941,
        5.878906683738906e-39,
        0.865003,
        0.865003,
        0.865003,
        1.1757813367477812e-38,
        0.705882,
        0.0156863,
        0.14902,
    ]
    substanceConcentrationLUT.UseLogScale = 0
    substanceConcentrationLUT.UseOpacityControlPointsFreehandDrawing = 0
    substanceConcentrationLUT.ShowDataHistogram = 0
    substanceConcentrationLUT.AutomaticDataHistogramComputation = 0
    substanceConcentrationLUT.DataHistogramNumberOfBins = 10
    substanceConcentrationLUT.ColorSpace = "Diverging"
    substanceConcentrationLUT.UseBelowRangeColor = 0
    substanceConcentrationLUT.BelowRangeColor = [0.0, 0.0, 0.0]
    substanceConcentrationLUT.UseAboveRangeColor = 0
    substanceConcentrationLUT.AboveRangeColor = [0.5, 0.5, 0.5]
    substanceConcentrationLUT.NanColor = [1.0, 1.0, 0.0]
    substanceConcentrationLUT.NanOpacity = 1.0
    substanceConcentrationLUT.Discretize = 1
    substanceConcentrationLUT.NumberOfTableValues = 256
    substanceConcentrationLUT.ScalarRangeInitialized = 1.0
    substanceConcentrationLUT.HSVWrap = 0
    substanceConcentrationLUT.VectorComponent = 0
    substanceConcentrationLUT.VectorMode = "Magnitude"
    substanceConcentrationLUT.AllowDuplicateScalars = 1
    substanceConcentrationLUT.Annotations = []
    substanceConcentrationLUT.ActiveAnnotatedValues = []
    substanceConcentrationLUT.IndexedColors = []
    substanceConcentrationLUT.IndexedOpacities = []

    # get opacity transfer function/opacity map for 'SubstanceConcentration'
    substanceConcentrationPWF = GetOpacityTransferFunction(
        "SubstanceConcentration"
    )
    substanceConcentrationPWF.Points = [
        0.0,
        0.0,
        0.5,
        0.0,
        1.1757813367477812e-38,
        1.0,
        0.5,
        0.0,
    ]
    substanceConcentrationPWF.AllowDuplicateScalars = 1
    substanceConcentrationPWF.UseLogScale = 0
    substanceConcentrationPWF.ScalarRangeInitialized = 1

    # Properties modified on vEGFconcentration
    vEGFconcentration.TimeArray = "None"

    # update the view to ensure updated data information
    renderView1.Update()

    # Rescale transfer function
    substanceConcentrationLUT.RescaleTransferFunction(0.0, 0.9964508973497085)

    # Rescale transfer function
    substanceConcentrationPWF.RescaleTransferFunction(0.0, 0.9964508973497085)

    # Apply a preset using its name. Note this may not work as expected when presets have duplicate names.
    substanceConcentrationLUT.ApplyPreset("Black, Blue and White", True)

    # invert the transfer function
    substanceConcentrationLUT.InvertTransferFunction()

    # Properties modified on substanceConcentrationPWF
    substanceConcentrationPWF.Points = [
        0.0,
        0.0,
        0.5,
        0.0,
        0.06396137177944183,
        0.0,
        0.5,
        0.0,
        0.08752609044313431,
        0.06951871514320374,
        0.5,
        0.0,
        0.11445719748735428,
        0.0,
        0.5,
        0.0,
        0.17168579995632172,
        0.0,
        0.5,
        0.0,
        0.21881522238254547,
        0.08021390438079834,
        0.5,
        0.0,
        0.2726774513721466,
        0.0,
        0.5,
        0.0,
        0.31307411193847656,
        0.0,
        0.5,
        0.0,
        0.35347074270248413,
        0.11229946464300156,
        0.5,
        0.0,
        0.3972338140010834,
        0.0,
        0.5,
        0.0,
        0.461195170879364,
        0.0,
        0.5,
        0.0,
        0.525156557559967,
        0.16042780876159668,
        0.5,
        0.0,
        0.5857515335083008,
        0.005347593687474728,
        0.5,
        0.0,
        0.6463465094566345,
        0.0,
        0.5,
        0.0,
        0.6901095509529114,
        0.1978609710931778,
        0.5,
        0.0,
        0.7406054139137268,
        0.0,
        0.5,
        0.0,
        0.8180323243141174,
        0.0,
        0.5,
        0.0,
        0.8550626039505005,
        0.32620322704315186,
        0.5,
        0.0,
        0.8887264728546143,
        0.0,
        0.5,
        0.0,
        0.9358559250831604,
        0.0,
        0.5,
        0.0,
        0.9964508973497085,
        0.614973247051239,
        0.5,
        0.0,
    ]

    # find source
    vessels = FindSource("Vessels")

    # set active source
    SetActiveSource(vessels)

    # get display properties
    vesselsDisplay = GetDisplayProperties(vessels, view=renderView1)
    renderView1.OrientationAxesVisibility = show_orientation_axes

    # change solid color
    vesselsDisplay.AmbientColor = [0.3333333333333333, 0.0, 0.0]
    vesselsDisplay.DiffuseColor = [0.3333333333333333, 0.0, 0.0]

    # get camera animation track for the view
    cameraAnimationCue1 = GetCameraTrack(view=renderView1)

    # create keyframes for this animation track

    # create a key frame
    keyFrame13766 = CameraKeyFrame()
    keyFrame13766.KeyTime = 0.0
    keyFrame13766.KeyValues = [0.0]
    keyFrame13766.Position = [
        954.9620344925504,
        2095.70409319902,
        371.5019581517663,
    ]
    keyFrame13766.FocalPoint = [0.0, 4.5, 4.5]
    keyFrame13766.ViewUp = [
        -0.16692216632147164,
        -0.0959876499942075,
        0.9812865847646836,
    ]
    keyFrame13766.ViewAngle = 30.0
    keyFrame13766.ParallelScale = 602.5414923472076
    keyFrame13766.PositionPathPoints = [
        954.962,
        2095.7,
        371.502,
        -1036.5033301689834,
        2088.8849682806476,
        32.076801886572014,
        -2247.4606194372627,
        512.4854799059449,
        -328.1136583656836,
        -1766.0343299588153,
        -1446.4374976886072,
        -437.8386176404253,
        45.251722982468245,
        -2312.775007627799,
        -214.4729436225022,
        1822.4621980197176,
        -1434.1572348627974,
        173.78445532319407,
        2227.321358962914,
        527.7987171780193,
        434.5678387829498,
    ]
    keyFrame13766.FocalPathPoints = [0.0, 4.5, 4.5]
    keyFrame13766.PositionMode = "Path"
    keyFrame13766.FocalPointMode = "Path"
    keyFrame13766.ClosedFocalPath = 0
    keyFrame13766.ClosedPositionPath = 1

    # create a key frame
    keyFrame13767 = CameraKeyFrame()
    keyFrame13767.KeyTime = 1.0
    keyFrame13767.KeyValues = [0.0]
    keyFrame13767.Position = [
        954.9620344925504,
        2095.70409319902,
        371.5019581517663,
    ]
    keyFrame13767.FocalPoint = [0.0, 4.5, 4.5]
    keyFrame13767.ViewUp = [
        -0.16692216632147164,
        -0.0959876499942075,
        0.9812865847646836,
    ]
    keyFrame13767.ViewAngle = 30.0
    keyFrame13767.ParallelScale = 602.5414923472076
    keyFrame13767.PositionPathPoints = [
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
    keyFrame13767.FocalPathPoints = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0]
    keyFrame13767.PositionMode = "Path"
    keyFrame13767.FocalPointMode = "Path"
    keyFrame13767.ClosedFocalPath = 0
    keyFrame13767.ClosedPositionPath = 0

    # initialize the animation track
    cameraAnimationCue1.TimeMode = "Normalized"
    cameraAnimationCue1.StartTime = 0.0
    cameraAnimationCue1.EndTime = 1.0
    cameraAnimationCue1.Enabled = 1
    cameraAnimationCue1.Mode = "Path-based"
    cameraAnimationCue1.Interpolation = "Spline"
    cameraAnimationCue1.KeyFrames = [keyFrame13766, keyFrame13767]
    cameraAnimationCue1.DataSource = None

    # get layout
    layout1 = GetLayout()

    # layout/tab size in pixels
    layout1.SetSize(952, 854)

    # current camera placement for renderView1
    renderView1.CameraPosition = [
        954.9619999999993,
        2095.7000000000003,
        371.50199999999995,
    ]
    renderView1.CameraFocalPoint = [0.0, 4.5, 4.5]
    renderView1.CameraViewUp = [
        -0.16692216632147164,
        -0.0959876499942075,
        0.9812865847646836,
    ]
    renderView1.CameraParallelScale = 602.5414923472076

    # get the absolute path of the filename
    abs = os.path.abspath(filename)
    # get the directory of the filename
    dir = os.path.dirname(abs)

    # save animation
    print("<pvpython> Create folder ..")
    output_folder_1 = os.path.join(dir, "vessel_vegf")
    output_folder_2 = os.path.join(dir, "vessel")
    output_folder_1 = output_folder_1 + "_bg{}_cb{}_ax{}".format(
        transparent_background, int(show_color_bar), show_orientation_axes
    )
    output_folder_2 = output_folder_2 + "_bg{}_cb{}_ax{}".format(
        transparent_background, int(show_color_bar), show_orientation_axes
    )
    output_folder_1 = os.path.abspath(output_folder_1)
    output_folder_2 = os.path.abspath(output_folder_2)

    if os.path.exists(output_folder_1):
        if not overwrite:
            print("<pvpython> folder 'vessel_vegf' already exists, aborting")
            return
        shutil.rmtree(output_folder_1)
    if os.path.exists(output_folder_2):
        if not overwrite:
            print("<pvpython> folder 'vessel' already exists, aborting")
            return
        shutil.rmtree(output_folder_2)
    os.makedirs(output_folder_1)
    os.makedirs(output_folder_2)

    # layout/tab size in pixels
    layout1.SetSize(952, 854)

    # current camera placement for renderView1
    renderView1.CameraPosition = [
        954.9619999999993,
        2095.7000000000003,
        371.50199999999995,
    ]
    renderView1.CameraFocalPoint = [0.0, 4.5, 4.5]
    renderView1.CameraViewUp = [
        -0.16692216632147164,
        -0.0959876499942075,
        0.9812865847646836,
    ]
    renderView1.CameraParallelScale = 602.5414923472076

    # save animation
    print("<pvpython> Save animation 1 ..")
    SaveAnimation(
        os.path.join(output_folder_1, "img.png"),
        renderView1,
        # ImageResolution=[952, 854],
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

    # hide data in view
    Hide(vEGFconcentration, renderView1)

    # layout/tab size in pixels
    layout1.SetSize(952, 854)

    # current camera placement for renderView1
    renderView1.CameraPosition = [
        954.9619999999993,
        2095.7000000000003,
        371.50199999999995,
    ]
    renderView1.CameraFocalPoint = [0.0, 4.5, 4.5]
    renderView1.CameraViewUp = [
        -0.16692216632147164,
        -0.0959876499942075,
        0.9812865847646836,
    ]
    renderView1.CameraParallelScale = 602.5414923472076

    # save animation
    print("<pvpython> Save animation 2 ..")
    SaveAnimation(
        os.path.join(output_folder_2, "img.png"),
        renderView1,
        # ImageResolution=[952, 854],
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
    # layout1.SetSize(952, 854)

    # # -----------------------------------
    # # saving camera placements for views

    # # current camera placement for renderView1
    # renderView1.CameraPosition = [-812.383, 174.638, 2179.52]
    # renderView1.CameraFocalPoint = [0.0, 4.5, 4.5]
    # renderView1.CameraViewUp = [
    #     0.01711090958076122,
    #     0.997296151690281,
    #     -0.07146749328943233,
    # ]
    # renderView1.CameraParallelScale = 602.5414923472076

    # # --------------------------------------------
    # # uncomment the following to render all views
    # # RenderAllViews()
    # # alternatively, if you want to write images, you can use SaveScreenshot(...).


def main(argc, argv):
    if argc != 6:
        print(
            "Usage: visualize.py <filename> <transparent_background> "
            + "<show_orientation_axes> <show_color_bar> <overwrite>"
        )
        return
    filename = argv[1]
    transparent_background = int(argv[2])
    show_orientation_axes = int(argv[3])
    show_color_bar = bool(int(argv[4]))
    overwrite = bool(int(argv[5]))
    # check if file filename exists
    if not os.path.isfile(filename):
        print("File {} does not exist".format(filename))
        return
    visualize(
        filename,
        transparent_background,
        show_orientation_axes,
        show_color_bar,
        overwrite,
    )


if __name__ == "__main__":
    main(len(sys.argv), sys.argv)
