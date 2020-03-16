# trace generated using paraview version 5.5.0

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

#-------------------------------------------------------------------------------------------------------------------------
# create a new 'XDMF Reader'
pfd_xdmf = XDMFReader(FileNames=['/disk/plasmaz/jaa/Critical_Balance/CB8waves_104_1/Raw_outputs/pfd.000400.xdmf']) #CHANGE THE FILE
#pfd_xdmf = XDMFReader(FileNames=['/disk/plasmaz/jaa/Critical_Balance/CB8waves_104_1/pfd.002000.xdmf'])

# Properties modified on pfd000400xdmf
pfd_xdmf.CellArrayStatus = ['hx', 'hy', 'hz', 'ex', 'ey', 'ez']
pfd_xdmf.GridStatus = ['patch-mrc_domain-0']

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [1550, 803]

# show data in view
pfd_xdmfDisplay = Show(pfd_xdmf, renderView1)

# trace defaults for the display properties.
pfd_xdmfDisplay.Representation = 'Outline'

# reset view to fit data
renderView1.ResetCamera()

# update the view to ensure updated data information
renderView1.Update()
#-------------------------------------------------------------------------------------------------------------------------


#-------------------------------------------------------------------------------------------------------------------------
# create a new 'Calculator'
calculator1 = Calculator(Input=pfd_xdmf)

# rename source object
RenameSource('B_field', calculator1)

# Properties modified on calculator1
calculator1.AttributeType = 'Cell Data'
calculator1.ResultArrayName = 'B_field'
calculator1.Function = 'iHat*hx+jHat*hy+kHat*hz'

# show data in view
calculator1Display = Show(calculator1, renderView1)

# trace defaults for the display properties.
calculator1Display.Representation = 'Outline'

# hide data in view
Hide(pfd_xdmf, renderView1)

# update the view to ensure updated data information
renderView1.Update()

# set active source
SetActiveSource(pfd_xdmf)

# create a new 'Calculator'
calculator1_1 = Calculator(Input=pfd_xdmf)

# Properties modified on calculator1_1
calculator1_1.AttributeType = 'Cell Data'
calculator1_1.ResultArrayName = 'Epar'
calculator1_1.Function = '(ex*hx + ey*hy + ez*hz)/sqrt(hx*hx +hy*hy +hz*hz)'

# show data in view
calculator1_1Display = Show(calculator1_1, renderView1)

# trace defaults for the display properties.
calculator1_1Display.Representation = 'Outline'

# hide data in view
Hide(pfd_xdmf, renderView1)

# update the view to ensure updated data information
renderView1.Update()

# rename source object
RenameSource('Epar', calculator1_1)
#-------------------------------------------------------------------------------------------------------------------------


#-------------------------------------------------------------------------------------------------------------------------
# set active source
SetActiveSource(calculator1)

# create a new 'Stream Tracer'
streamTracer1 = StreamTracer(Input=calculator1,
    SeedType='High Resolution Line Source')

# toggle 3D widget visibility (only when running from the GUI)
Hide3DWidgets(proxy=streamTracer1.SeedType)

# Properties modified on streamTracer1.SeedType
streamTracer1.SeedType.Point1 = [6.0, 6.0, 0.0]
streamTracer1.SeedType.Point2 = [6.0, 6.0, 125.0000638961792]
streamTracer1.SeedType.Resolution = 50

# Properties modified on streamTracer1
streamTracer1.IntegratorType = 'Runge-Kutta 4'

# Properties modified on streamTracer1.SeedType
streamTracer1.SeedType.Point1 = [6.0, 6.0, 0.0]
streamTracer1.SeedType.Point2 = [6.0, 6.0, 125.0000638961792]
streamTracer1.SeedType.Resolution = 50

# show data in view
streamTracer1Display = Show(streamTracer1, renderView1)

# trace defaults for the display properties.
streamTracer1Display.Representation = 'Surface'

# update the view to ensure updated data information
renderView1.Update()

# change representation type
streamTracer1Display.SetRepresentationType('3D Glyphs')

# set scalar coloring
ColorBy(streamTracer1Display, ('POINTS', 'B_field', 'Magnitude'))

# rescale color and/or opacity maps used to include current data range
streamTracer1Display.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
streamTracer1Display.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'B_field'
b_fieldLUT = GetColorTransferFunction('B_field')

# Properties modified on streamTracer1Display
streamTracer1Display.Scaling = 1

# Properties modified on streamTracer1Display
streamTracer1Display.ScaleMode = 'Magnitude'

# Properties modified on streamTracer1Display
streamTracer1Display.ScaleFactor = 1.0

# Properties modified on streamTracer1Display
streamTracer1Display.SelectScaleArray = 'B_field'

# Properties modified on streamTracer1Display
streamTracer1Display.GlyphType = 'Sphere'

# set active source
SetActiveSource(calculator1)

# create a new 'Stream Tracer'
streamTracer2 = StreamTracer(Input=calculator1,
    SeedType='High Resolution Line Source')

# toggle 3D widget visibility (only when running from the GUI)
Hide3DWidgets(proxy=streamTracer2.SeedType)

# Properties modified on streamTracer2.SeedType
streamTracer2.SeedType.Point1 = [6.0, 18.0, 0.0]
streamTracer2.SeedType.Point2 = [6.0, 18.0, 125.0000638961792]
streamTracer2.SeedType.Resolution = 50

# Properties modified on streamTracer2
streamTracer2.IntegratorType = 'Runge-Kutta 4'

# Properties modified on streamTracer2.SeedType
streamTracer2.SeedType.Point1 = [6.0, 18.0, 0.0]
streamTracer2.SeedType.Point2 = [6.0, 18.0, 125.0000638961792]
streamTracer2.SeedType.Resolution = 50

# show data in view
streamTracer2Display = Show(streamTracer2, renderView1)

# trace defaults for the display properties.
streamTracer2Display.Representation = 'Surface'

# update the view to ensure updated data information
renderView1.Update()

# change representation type
streamTracer2Display.SetRepresentationType('3D Glyphs')

# set scalar coloring
ColorBy(streamTracer2Display, ('POINTS', 'B_field', 'Magnitude'))

# rescale color and/or opacity maps used to include current data range
streamTracer2Display.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
streamTracer2Display.SetScalarBarVisibility(renderView1, True)

# Properties modified on streamTracer2Display
streamTracer2Display.Scaling = 1

# Properties modified on streamTracer2Display
streamTracer2Display.ScaleMode = 'Magnitude'

# Properties modified on streamTracer2Display
streamTracer2Display.ScaleFactor = 1.0

# Properties modified on streamTracer2Display
streamTracer2Display.SelectScaleArray = 'B_field'

# Properties modified on streamTracer2Display
streamTracer2Display.GlyphType = 'Sphere'

# set active source
SetActiveSource(calculator1)

# create a new 'Stream Tracer'
streamTracer3 = StreamTracer(Input=calculator1,
    SeedType='High Resolution Line Source')

# toggle 3D widget visibility (only when running from the GUI)
Hide3DWidgets(proxy=streamTracer3.SeedType)

# toggle 3D widget visibility (only when running from the GUI)
Show3DWidgets(proxy=streamTracer3.SeedType)

# toggle 3D widget visibility (only when running from the GUI)
Hide3DWidgets(proxy=streamTracer3.SeedType)

# toggle 3D widget visibility (only when running from the GUI)
Show3DWidgets(proxy=streamTracer3.SeedType)

# toggle 3D widget visibility (only when running from the GUI)
Hide3DWidgets(proxy=streamTracer3.SeedType)

# Properties modified on streamTracer3.SeedType
streamTracer3.SeedType.Point1 = [18.0, 6.0, 0.0]
streamTracer3.SeedType.Point2 = [18.0, 6.0, 125.0000638961792]
streamTracer3.SeedType.Resolution = 50

# Properties modified on streamTracer3
streamTracer3.IntegratorType = 'Runge-Kutta 4'

# Properties modified on streamTracer3.SeedType
streamTracer3.SeedType.Point1 = [18.0, 6.0, 0.0]
streamTracer3.SeedType.Point2 = [18.0, 6.0, 125.0000638961792]
streamTracer3.SeedType.Resolution = 50

# show data in view
streamTracer3Display = Show(streamTracer3, renderView1)

# trace defaults for the display properties.
streamTracer3Display.Representation = 'Surface'

# update the view to ensure updated data information
renderView1.Update()

# change representation type
streamTracer3Display.SetRepresentationType('3D Glyphs')

# set scalar coloring
ColorBy(streamTracer3Display, ('POINTS', 'B_field', 'Magnitude'))

# rescale color and/or opacity maps used to include current data range
streamTracer3Display.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
streamTracer3Display.SetScalarBarVisibility(renderView1, True)

# Properties modified on streamTracer3Display
streamTracer3Display.Scaling = 1

# Properties modified on streamTracer3Display
streamTracer3Display.ScaleMode = 'Magnitude'

# Properties modified on streamTracer3Display
streamTracer3Display.ScaleFactor = 1.0

# Properties modified on streamTracer3Display
streamTracer3Display.SelectScaleArray = 'B_field'

# Properties modified on streamTracer3Display
streamTracer3Display.GlyphType = 'Sphere'

# set active source
SetActiveSource(calculator1)

# create a new 'Stream Tracer'
streamTracer4 = StreamTracer(Input=calculator1,
    SeedType='High Resolution Line Source')

# toggle 3D widget visibility (only when running from the GUI)
Hide3DWidgets(proxy=streamTracer4.SeedType)

# Properties modified on streamTracer4.SeedType
streamTracer4.SeedType.Point1 = [18.0, 18.0, 0.0]
streamTracer4.SeedType.Point2 = [18.0, 18.0, 125.0000638961792]
streamTracer4.SeedType.Resolution = 50

# Properties modified on streamTracer4
streamTracer4.IntegratorType = 'Runge-Kutta 4'

# Properties modified on streamTracer4.SeedType
streamTracer4.SeedType.Point1 = [18.0, 18.0, 0.0]
streamTracer4.SeedType.Point2 = [18.0, 18.0, 125.0000638961792]
streamTracer4.SeedType.Resolution = 50

# show data in view
streamTracer4Display = Show(streamTracer4, renderView1)

# trace defaults for the display properties.
streamTracer4Display.Representation = 'Surface'

# update the view to ensure updated data information
renderView1.Update()

# change representation type
streamTracer4Display.SetRepresentationType('3D Glyphs')

# set scalar coloring
ColorBy(streamTracer4Display, ('POINTS', 'B_field', 'Magnitude'))

# rescale color and/or opacity maps used to include current data range
streamTracer4Display.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
streamTracer4Display.SetScalarBarVisibility(renderView1, True)

# Properties modified on streamTracer4Display
streamTracer4Display.Scaling = 1

# Properties modified on streamTracer4Display
streamTracer4Display.ScaleMode = 'Magnitude'

# Properties modified on streamTracer4Display
streamTracer4Display.ScaleFactor = 1.0

# Properties modified on streamTracer4Display
streamTracer4Display.SelectScaleArray = 'B_field'

# Properties modified on streamTracer4Display
streamTracer4Display.GlyphType = 'Sphere'

#-------------------------------------------------------------------------------------------------------------------------


#-------------------------------------------------------------------------------------------------------------------------

# set active source
SetActiveSource(calculator1_1)

# set scalar coloring
ColorBy(calculator1_1Display, ('CELLS', 'Epar'))

# rescale color and/or opacity maps used to include current data range
calculator1_1Display.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
calculator1_1Display.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'Epar'
eparLUT = GetColorTransferFunction('Epar')

# change representation type
calculator1_1Display.SetRepresentationType('Volume')

# Properties modified on calculator1_1Display
calculator1_1Display.ShowIsosurfaces = 1

calculator1_1Display.IsosurfaceValues = [-0.0.01441, 0.01441]
#-------------------------------------------------------------------------------------------------------------------------


#-------------------------------------------------------------------------------------------------------------------------

# get color transfer function/color map for 'Epar'
eparLUT = GetColorTransferFunction('Epar')

# Rescale transfer function
eparLUT.RescaleTransferFunction(-0.015, 0.015) #(0.001, 0.8)

# get opacity transfer function/opacity map for 'Epar'
eparPWF = GetOpacityTransferFunction('Epar')

# Rescale transfer function
eparPWF.RescaleTransferFunction(-0.015, 0.015) # (0.001, 0.8)

# find source
streamTracer1 = FindSource('StreamTracer1')

# set active source
SetActiveSource(streamTracer1)

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [1212, 528]

# get display properties
streamTracer1Display = GetDisplayProperties(streamTracer1, view=renderView1)

# get color transfer function/color map for 'B_field'
b_fieldLUT = GetColorTransferFunction('B_field')

# Rescale transfer function
b_fieldLUT.RescaleTransferFunction(0.015, 0.2)

# get opacity transfer function/opacity map for 'B_field'
b_fieldPWF = GetOpacityTransferFunction('B_field')

# Rescale transfer function
b_fieldPWF.RescaleTransferFunction(0.015, 0.2)

#-------------------------------------------------------------------------------------------------------------------------



#-------------------------------------------------------------------------------------------------------------------------
# Apply a preset using its name. Note this may not work as expected when presets have duplicate names.
eparLUT.ApplyPreset('Cool to Warm (Extended)', True)

# set active source
SetActiveSource(streamTracer1)

# Apply a preset using its name. Note this may not work as expected when presets have duplicate names.
b_fieldLUT.ApplyPreset('Black-Body Radiation', True)

# hide color bar/color legend
streamTracer1Display.SetScalarBarVisibility(renderView1, False)

# set active source
SetActiveSource(calculator1_1)

# hide color bar/color legend
calculator1_1Display.SetScalarBarVisibility(renderView1, False)

# Properties modified on renderView1
renderView1.UseGradientBackground = 1

# Properties modified on renderView1
renderView1.Background = [1.0, 1.0, 1.0]

# reset view to fit data
renderView1.ResetCamera()

# current camera placement for renderView1
renderView1.CameraPosition = [-129.59156150296562, 12.00055973418057, 62.50025046057999]
renderView1.CameraFocalPoint = [11.998267595332418, 12.00055973418057, 62.50025046057999]
renderView1.CameraViewUp = [0.0, -1.0, 2.220446049250313e-16]
renderView1.CameraParallelScale = 64.92088015465893

# save screenshot
SaveScreenshot('/disk/plasmaz/jaa/Critical_Balance/CB8waves_104_1/Video_maglines/Blines_Epar_mzmy_1.png', renderView1, ImageResolution=[1550, 803]) #HERE SAVE THE PIC

# reset view to fit data
renderView1.ResetCamera()

# current camera placement for renderView1
renderView1.CameraPosition = [11.998267595332418, 153.5903888324786, 62.50025046057999]
renderView1.CameraFocalPoint = [11.998267595332418, 12.00055973418057, 62.50025046057999]
renderView1.CameraViewUp = [-1.0, 0.0, 2.220446049250313e-16]
renderView1.CameraParallelScale = 64.92088015465893

# save screenshot
SaveScreenshot('/disk/plasmaz/jaa/Critical_Balance/CB8waves_104_1/Video_maglines/Blines_Epar_mzmx_1.png', renderView1, ImageResolution=[1550, 803]) #HERE SAVE THE PIC

# reset view to fit data
renderView1.ResetCamera()

# current camera placement for renderView1
renderView1.CameraPosition = [11.998267595332418, 12.00055973418057, 179.5166381451238]
renderView1.CameraFocalPoint = [11.998267595332418, 12.00055973418057, 62.50025046057999]
renderView1.CameraViewUp = [4.440892098500626e-16, -1.0, 0.0]
renderView1.CameraParallelScale = 64.92088015465893

# save screenshot
SaveScreenshot('/disk/plasmaz/jaa/Critical_Balance/CB8waves_104_1/Video_maglines/Blines_Epar_mxmy_1.png', renderView1, ImageResolution=[1550, 803]) #HERE SAVE THE PIC

#### saving camera placements for all active views

# current camera placement for renderView1
renderView1.CameraPosition = [11.998267595332418, 12.00055973418057, 179.5166381451238]
renderView1.CameraFocalPoint = [11.998267595332418, 12.00055973418057, 62.50025046057999]
renderView1.CameraViewUp = [4.440892098500626e-16, -1.0, 0.0]
renderView1.CameraParallelScale = 64.92088015465893

#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).