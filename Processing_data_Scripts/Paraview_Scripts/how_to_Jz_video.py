# trace generated using paraview version 5.5.0

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

#-------------------------------------------------------------------------------------------------------------------------
# create a new 'XDMF Reader'
pfd_xdmf = XDMFReader(FileNames=['/disk/plasmaz/jaa/Critical_Balance/CB8waves_104_1/Raw_outputs/pfd.010000.xdmf']) #CHANGE THE FILE
#pfd_xdmf = XDMFReader(FileNames=['/disk/plasmaz/jaa/Critical_Balance/CB8waves_104_1/pfd.002000.xdmf'])

# Properties modified on pfd_xdmf
pfd_xdmf.CellArrayStatus = ['jz']
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

# set scalar coloring
ColorBy(pfd_xdmfDisplay, ('CELLS', 'jz'))
# rescale color and/or opacity maps used to include current data range
pfd_xdmfDisplay.RescaleTransferFunctionToDataRange(True, False)
# show color bar/color legend
pfd_xdmfDisplay.SetScalarBarVisibility(renderView1, True)
# get color transfer function/color map for 'jz'
jzLUT = GetColorTransferFunction('jz')
# change representation type
pfd_xdmfDisplay.SetRepresentationType('Volume')

# Apply a preset using its name. Note this may not work as expected when presets have duplicate names.
jzLUT.ApplyPreset('Cool to Warm (Extended)', True)
# get color transfer function/color map for 'Jmag'
#jmagLUT.RescaleTransferFunction(0.2, 0.6) #(0.001, 0.8)
jzLUT.RescaleTransferFunction(-0.7, 0.7) #(0.001, 0.8)
# get opacity transfer function/opacity map for 'Jmag'
jzPWF = GetOpacityTransferFunction('Jmag')
# Rescale transfer function
jzPWF.RescaleTransferFunction(-0.7, 0.7) # (0.001, 0.8)
#------------------------------------------------------------------------------------------------------------------------

#jmagLUT = GetColorTransferFunction('Jmag')
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
SaveScreenshot('/disk/plasmaz/jaa/Critical_Balance/CB8waves_104_1/Video_maglines/Jz_mzmy_25.png', renderView1, ImageResolution=[1550, 803]) #HERE SAVE THE PIC

# reset view to fit data
renderView1.ResetCamera()

# current camera placement for renderView1
renderView1.CameraPosition = [11.998267595332418, 153.5903888324786, 62.50025046057999]
renderView1.CameraFocalPoint = [11.998267595332418, 12.00055973418057, 62.50025046057999]
renderView1.CameraViewUp = [-1.0, 0.0, 2.220446049250313e-16]
renderView1.CameraParallelScale = 64.92088015465893

# save screenshot
SaveScreenshot('/disk/plasmaz/jaa/Critical_Balance/CB8waves_104_1/Video_maglines/Jz_mzmx_25.png', renderView1, ImageResolution=[1550, 803]) #HERE SAVE THE PIC

# reset view to fit data
renderView1.ResetCamera()

# current camera placement for renderView1
renderView1.CameraPosition = [11.998267595332418, 12.00055973418057, 179.5166381451238]
renderView1.CameraFocalPoint = [11.998267595332418, 12.00055973418057, 62.50025046057999]
renderView1.CameraViewUp = [4.440892098500626e-16, -1.0, 0.0]
renderView1.CameraParallelScale = 64.92088015465893

# save screenshot
SaveScreenshot('/disk/plasmaz/jaa/Critical_Balance/CB8waves_104_1/Video_maglines/Jz_mxmy_25.png', renderView1, ImageResolution=[1550, 803]) #HERE SAVE THE PIC

#### saving camera placements for all active views

# current camera placement for renderView1
renderView1.CameraPosition = [11.998267595332418, 12.00055973418057, 179.5166381451238]
renderView1.CameraFocalPoint = [11.998267595332418, 12.00055973418057, 62.50025046057999]
renderView1.CameraViewUp = [4.440892098500626e-16, -1.0, 0.0]
renderView1.CameraParallelScale = 64.92088015465893



#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).