# trace generated using paraview version 5.5.0

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

#-------------------------------------------------------------------------------------------------------------------------
# create a new 'XDMF Reader'
#pfd_xdmf = XDMFReader(FileNames=['/disk/plasmaz/jaa/Critical_Balance/CB8waves_104_1/Raw_outputs/pfd.002000.xdmf']) #CHANGE THE FILE
pfd_xdmf = XDMFReader(FileNames=['/disk/plasmaz/jaa/Critical_Balance/CB8waves_104_1/pfd.002000.xdmf'])

# Properties modified on pfd000400xdmf
pfd_xdmf.CellArrayStatus = ['ex', 'ey', 'ez','hx', 'hy', 'hz', 'jx', 'jy', 'jz']
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
RenameSource('Bmag', calculator1)
# Properties modified on calculator1
calculator1.AttributeType = 'Cell Data'
calculator1.ResultArrayName = 'Bmag'
calculator1.Function = 'sqrt(hx*hx +hy*hy +hz*hz)'
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
calculator2 = Calculator(Input=pfd_xdmf)
# Properties modified on calculator2
calculator2.AttributeType = 'Cell Data'
calculator2.ResultArrayName = 'Jmag'
calculator2.Function = 'sqrt(jx*jx +jy*jy +jz*jz)'
# show data in view
calculator2Display = Show(calculator2, renderView1)
# trace defaults for the display properties.
calculator2Display.Representation = 'Outline'
# hide data in view
Hide(pfd_xdmf, renderView1)
# update the view to ensure updated data information
renderView1.Update()
# rename source object
RenameSource('Jmag', calculator2)


# set active source
SetActiveSource(pfd_xdmf)
# create a new 'Calculator'
calculator3 = Calculator(Input=pfd_xdmf)
# Properties modified on calculator2
calculator3.AttributeType = 'Cell Data'
calculator3.ResultArrayName = 'Epar'
calculator3.Function = '(ex*hx+ey*hy+ez*hz)/(sqrt(hx*hx+hy*hy+hz*hz))'
# show data in view
calculator3Display = Show(calculator3, renderView1)
# trace defaults for the display properties.
calculator3Display.Representation = 'Outline'
# hide data in view
Hide(pfd_xdmf, renderView1)
# update the view to ensure updated data information
renderView1.Update()
# rename source object
RenameSource('Epar', calculator3)


SetActiveSource(pfd_xdmf)
# create a new 'Calculator'
calculator1_3 = Calculator(Input=pfd_xdmf)
# Properties modified on calculator2
calculator1_3.Function = ''
# show data in view
calculator1_3Display = Show(calculator1_3, renderView1)
# trace defaults for the display properties.
calculator1_3Display.Representation = 'Outline'

#-------------------------------------------------------------------------------------------------------------------------

#B
#-------------------------------------------------------------------------------------------------------------------------
# set active source
SetActiveSource(calculator1)
# set scalar coloring
ColorBy(calculator1Display, ('CELLS', 'Bmag'))
# rescale color and/or opacity maps used to include current data range
calculator1Display.RescaleTransferFunctionToDataRange(True, False)
# show color bar/color legend
calculator1Display.SetScalarBarVisibility(renderView1, True)
# get color transfer function/color map for 'Jmag'
jmagLUT = GetColorTransferFunction('Bmag')
# change representation type
calculator1Display.SetRepresentationType('Volume')
# Properties modified on calculator2Display
calculator1Display.ShowIsosurfaces = 1
calculator1Display.IsosurfaceValues = [0.02]
#-------------------------------------------------------------------------------------------------------------------------
#J
#-------------------------------------------------------------------------------------------------------------------------
# set active source
SetActiveSource(calculator2)
# set scalar coloring
ColorBy(calculator2Display, ('CELLS', 'Jmag'))
# rescale color and/or opacity maps used to include current data range
calculator2Display.RescaleTransferFunctionToDataRange(True, False)
# show color bar/color legend
calculator2Display.SetScalarBarVisibility(renderView1, True)
# get color transfer function/color map for 'Jmag'
jmagLUT = GetColorTransferFunction('Jmag')
# change representation type
calculator2Display.SetRepresentationType('Volume')
# Properties modified on calculator2Display
calculator2Display.ShowIsosurfaces = 1
calculator2Display.IsosurfaceValues = [0.5]
#-------------------------------------------------------------------------------------------------------------------------
#Epar
#-------------------------------------------------------------------------------------------------------------------------
# set active source
SetActiveSource(calculator3)
# set scalar coloring
ColorBy(calculator3Display, ('CELLS', 'Epar'))
# rescale color and/or opacity maps used to include current data range
calculator3Display.RescaleTransferFunctionToDataRange(True, False)
# show color bar/color legend
calculator3Display.SetScalarBarVisibility(renderView1, True)
# get color transfer function/color map for 'Epar'
eparLUT = GetColorTransferFunction('Epar')
# change representation type
calculator3Display.SetRepresentationType('Volume')
# Properties modified on calculator2Display
calculator3Display.ShowIsosurfaces = 1
calculator3Display.IsosurfaceValues = [0.01, -0.01]
#-------------------------------------------------------------------------------------------------------------------------



#-------------------------------------------------------------------------------------------------------------------------
# set active source
SetActiveSource(calculator2)
# get color transfer function/color map for 'Jmag'
jmagLUT = GetColorTransferFunction('Jmag')
# Rescale transfer function
jmagLUT.RescaleTransferFunction(0.001, 0.8)
# get opacity transfer function/opacity map for 'Jmag'
jmagPWF = GetOpacityTransferFunction('Jmag')
# Rescale transfer function
jmagPWF.RescaleTransferFunction(0.001, 0.8)
# Apply a preset using its name. Note this may not work as expected when presets have duplicate names.
jmagLUT.ApplyPreset('Black, Blue and White', True)


# set active source
SetActiveSource(calculator1)
# get active view
# get color transfer function/color map for 'Bmag'
bmagLUT = GetColorTransferFunction('Bmag')
# Rescale transfer function
bmagLUT.RescaleTransferFunction(0.05, 0.2)
# get opacity transfer function/opacity map for 'Bmag'
bmagPWF = GetOpacityTransferFunction('Bmag')
# Rescale transfer function
bmagPWF.RescaleTransferFunction(0.05, 0.2)
# Apply a preset using its name. Note this may not work as expected when presets have duplicate names.
bmagLUT.ApplyPreset('magenta', True)

# set active source
SetActiveSource(calculator3)
# get color transfer function/color map for 'Epar'
eparLUT = GetColorTransferFunction('Epar')
# Rescale transfer function
eparLUT.RescaleTransferFunction(0.001, 0.8)
# get opacity transfer function/opacity map for 'Epar'
eparPWF = GetOpacityTransferFunction('Epar')
# Rescale transfer function
eparPWF.RescaleTransferFunction(0.001, 0.8)
# Apply a preset using its name. Note this may not work as expected when presets have duplicate names.
eparLUT.ApplyPreset('Cool to Warm (Extended)', True)

#-------------------------------------------------------------------------------------------------------------------------


#-------------------------------------------------------------------------------------------------------------------------


# set active source
SetActiveSource(calculator1)

# hide color bar/color legend
calculator1Display.SetScalarBarVisibility(renderView1, False)

# set active source
SetActiveSource(calculator2)

# hide color bar/color legend
calculator2Display.SetScalarBarVisibility(renderView1, False)

# set active source
SetActiveSource(calculator3)

# hide color bar/color legend
calculator3Display.SetScalarBarVisibility(renderView1, False)


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
SaveScreenshot('/disk/plasmaz/jaa/Critical_Balance/CB8waves_104_1/Video_maglines/BEparJ_mzmy_5.png', renderView1, ImageResolution=[1550, 803]) #HERE SAVE THE PIC

# reset view to fit data
renderView1.ResetCamera()

# current camera placement for renderView1
renderView1.CameraPosition = [11.998267595332418, 153.5903888324786, 62.50025046057999]
renderView1.CameraFocalPoint = [11.998267595332418, 12.00055973418057, 62.50025046057999]
renderView1.CameraViewUp = [-1.0, 0.0, 2.220446049250313e-16]
renderView1.CameraParallelScale = 64.92088015465893

# save screenshot
SaveScreenshot('/disk/plasmaz/jaa/Critical_Balance/CB8waves_104_1/Video_maglines/BEparJ_mzmx_5.png', renderView1, ImageResolution=[1550, 803]) #HERE SAVE THE PIC

# reset view to fit data
renderView1.ResetCamera()

# current camera placement for renderView1
renderView1.CameraPosition = [11.998267595332418, 12.00055973418057, 179.5166381451238]
renderView1.CameraFocalPoint = [11.998267595332418, 12.00055973418057, 62.50025046057999]
renderView1.CameraViewUp = [4.440892098500626e-16, -1.0, 0.0]
renderView1.CameraParallelScale = 64.92088015465893

# save screenshot
SaveScreenshot('/disk/plasmaz/jaa/Critical_Balance/CB8waves_104_1/Video_maglines/BEparJ_mxmy_5.png', renderView1, ImageResolution=[1550, 803]) #HERE SAVE THE PIC

#### saving camera placements for all active views

# current camera placement for renderView1
renderView1.CameraPosition = [11.998267595332418, 12.00055973418057, 179.5166381451238]
renderView1.CameraFocalPoint = [11.998267595332418, 12.00055973418057, 62.50025046057999]
renderView1.CameraViewUp = [4.440892098500626e-16, -1.0, 0.0]
renderView1.CameraParallelScale = 64.92088015465893

#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).