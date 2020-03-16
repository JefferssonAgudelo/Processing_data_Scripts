# trace generated using paraview version 5.5.0

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

#-------------------------------------------------------------------------------------------------------------------------
# create a new 'XDMF Reader'
pfd_xdmf = XDMFReader(FileNames=['/disk/plasmaz/jaa/Critical_Balance/CB8waves_104_1/Raw_outputs/pfd.003600.xdmf']) #CHANGE THE FILE
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

#---------------------------------------------------------------

# set active source
SetActiveSource(pfd_xdmf)
# create a new 'Calculator'
calculator1_1 = Calculator(Input=pfd_xdmf)
# Properties modified on calculator1_1
calculator1_1.AttributeType = 'Cell Data'
calculator1_1.ResultArrayName = 'Eparmag'
calculator1_1.Function = '(ex*hx +ey*hy +ez*hz)/sqrt(hx*hx + hy*hy + hz*hz)'
# show data in view
calculator1_1Display = Show(calculator1_1, renderView1)
# trace defaults for the display properties.
calculator1_1Display.Representation = 'Outline'
# hide data in view
Hide(pfd_xdmf, renderView1)
# update the view to ensure updated data information
renderView1.Update()
# rename source object
RenameSource('Eparmag', calculator1_1)
#-------------------------------------------------------------------------------------------------------------------------


#-------------------------------------------------------------------------------------------------------------------------
# erdc_green2yellow_BW
# erdc_gold_BW
# RED_TEMPERATURE

# set active source
SetActiveSource(calculator1_1)
# set scalar coloring
ColorBy(calculator1_1Display, ('CELLS', 'Eparmag'))
# rescale color and/or opacity maps used to include current data range
calculator1_1Display.RescaleTransferFunctionToDataRange(True, False)
# show color bar/color legend
calculator1_1Display.SetScalarBarVisibility(renderView1, True)
# get color transfer function/color map for 'Eparmag'
eparLUT = GetColorTransferFunction('Eparmag')
# change representation type
calculator1_1Display.SetRepresentationType('Volume')
# Properties modified on calculator1_1Display
calculator1_1Display.ShowIsosurfaces = 1
calculator1_1Display.IsosurfaceValues = [-0.01046, 0.01046]
# hide color bar/color legend
calculator1_1Display.SetScalarBarVisibility(renderView1, False)
# Apply a preset using its name. Note this may not work as expected when presets have duplicate names.
eparLUT.ApplyPreset('Cool to Warm (Extended)', True)
# get color transfer function/color map for 'Eparmag'
eparLUT.RescaleTransferFunction(-0.015, 0.015) #(0.001, 0.8)
# get opacity transfer function/opacity map for 'Eparmag'
eparPWF = GetOpacityTransferFunction('Eparmag')
# Rescale transfer function
eparPWF.RescaleTransferFunction(-0.015, 0.015) # (0.001, 0.8)

#eparLUT = GetColorTransferFunction('Eparmag')
# Rescale transfer function
#-------------------------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------------------------------
# Properties modified on renderView1
renderView1.UseGradientBackground = 1
# Properties modified on renderView1
renderView1.Background = [1.0, 1.0, 1.0]

SetActiveSource(pfd_xdmf)
pfd_xdmfDisplay = Show(pfd_xdmf, renderView1)
#-------------------------------------------------------------------------------------------------------------------------


#-------------------------------------------------------------------------------------------------------------------------

# current camera placement for renderView1
renderView1.CameraPosition = [-129.59156150296562, 12.00055973418057, 62.50025046057999]
renderView1.CameraFocalPoint = [11.998267595332418, 12.00055973418057, 62.50025046057999]
renderView1.CameraViewUp = [0.0, -1.0, 2.220446049250313e-16]
renderView1.CameraParallelScale = 64.92088015465893

# save screenshot
SaveScreenshot('/disk/plasmaz/jaa/Critical_Balance/CB8waves_104_1/Video_maglines/Epar3_mzmy_9.png', renderView1, ImageResolution=[1550, 803]) #HERE SAVE THE PIC

# reset view to fit data
renderView1.ResetCamera()

# current camera placement for renderView1
renderView1.CameraPosition = [11.998267595332418, 153.5903888324786, 62.50025046057999]
renderView1.CameraFocalPoint = [11.998267595332418, 12.00055973418057, 62.50025046057999]
renderView1.CameraViewUp = [-1.0, 0.0, 2.220446049250313e-16]
renderView1.CameraParallelScale = 64.92088015465893

# save screenshot
SaveScreenshot('/disk/plasmaz/jaa/Critical_Balance/CB8waves_104_1/Video_maglines/Epar3_mzmx_9.png', renderView1, ImageResolution=[1550, 803]) #HERE SAVE THE PIC

# reset view to fit data
renderView1.ResetCamera()

# current camera placement for renderView1
renderView1.CameraPosition = [11.998267595332418, 12.00055973418057, 179.5166381451238]
renderView1.CameraFocalPoint = [11.998267595332418, 12.00055973418057, 62.50025046057999]
renderView1.CameraViewUp = [4.440892098500626e-16, -1.0, 0.0]
renderView1.CameraParallelScale = 64.92088015465893

# save screenshot
SaveScreenshot('/disk/plasmaz/jaa/Critical_Balance/CB8waves_104_1/Video_maglines/Epar3_mxmy_9.png', renderView1, ImageResolution=[1550, 803]) #HERE SAVE THE PIC

#### saving camera placements for all active views

# current camera placement for renderView1
renderView1.CameraPosition = [11.998267595332418, 12.00055973418057, 179.5166381451238]
renderView1.CameraFocalPoint = [11.998267595332418, 12.00055973418057, 62.50025046057999]
renderView1.CameraViewUp = [4.440892098500626e-16, -1.0, 0.0]
renderView1.CameraParallelScale = 64.92088015465893

#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).