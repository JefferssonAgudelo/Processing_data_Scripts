# trace generated using paraview version 5.5.0

#### import the simple module from the paraview
from paraview.simple import *
import glob, os
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

#path = '/disk/plasmaz/jaa/Critical_Balance/CB8waves_104_1/Raw_outputs'
path = '/disk/plasmaz/jaa/Critical_Balance/Temperature_test/This_are_the_runs_2/one'
file=sorted(glob.glob(os.path.join(path, '*.xdmf')))
a=-1

for i in glob.glob( os.path.join(path, '*.xdmf')):
  a=a+1;      
  filename1= file[a]#'/disk/plasma2/jaa/CB8WAVES/CB8waves_04/pfd.002000_p000000.h5' 
  #Read the file
  pfd_xdmf = XDMFReader(FileNames=filename1)
  # Properties modified on pfd000400xdmf
  pfd_xdmf.CellArrayStatus = ['hx', 'hy', 'hz', 'jx', 'jy', 'jz']
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
  calculator1.Function = 'sqrt(hx*hx + hy*hy + hz*hz)'
  # show data in view
  calculator1Display = Show(calculator1, renderView1)
  # trace defaults for the display properties.
  calculator1Display.Representation = 'Outline'
  # hide data in view
  Hide(pfd_xdmf, renderView1)
  # update the view to ensure updated data information
  renderView1.Update()
  
  #---------------------------------------------------------------
  
  # set active source
  SetActiveSource(pfd_xdmf)
  # create a new 'Calculator'
  calculator1_1 = Calculator(Input=pfd_xdmf)
  # Properties modified on calculator1_1
  calculator1_1.AttributeType = 'Cell Data'
  calculator1_1.ResultArrayName = 'Jmag'
  calculator1_1.Function = 'sqrt(jx*jx +jy*jy +jz*jz)'
  # show data in view
  calculator1_1Display = Show(calculator1_1, renderView1)
  # trace defaults for the display properties.
  calculator1_1Display.Representation = 'Outline'
  # hide data in view
  Hide(pfd_xdmf, renderView1)
  # update the view to ensure updated data information
  renderView1.Update()
  # rename source object
  RenameSource('Jmag', calculator1_1)
  #-------------------------------------------------------------------------------------------------------------------------
  
  
  #-------------------------------------------------------------------------------------------------------------------------
  # set active source
  SetActiveSource(calculator1)
  # set scalar coloring
  ColorBy(calculator1Display, ('CELLS', 'Bmag'))
  # rescale color and/or opacity maps used to include current data range
  calculator1Display.RescaleTransferFunctionToDataRange(True, False)
  # show color bar/color legend
  calculator1Display.SetScalarBarVisibility(renderView1, True)
  # get color transfer function/color map for 'Bmag'
  bmagLUT = GetColorTransferFunction('Bmag')
  # change representation type
  calculator1Display.SetRepresentationType('Volume')
  # Properties modified on calculator1Display
  calculator1Display.ShowIsosurfaces = 1
  # Properties modified on calculator1Display
  calculator1Display.IsosurfaceValues = [0.02]
  # hide color bar/color legend
  calculator1Display.SetScalarBarVisibility(renderView1, False)
  # Apply a preset using its name. Note this may not work as expected when presets have duplicate names.
  bmagLUT.ApplyPreset('Black-Body Radiation', True)
  # Rescale transfer function
  bmagLUT.RescaleTransferFunction(0.015, 0.025)
  # get opacity transfer function/opacity map for 'Bmag'
  bmagPWF = GetOpacityTransferFunction('Bmag')
  # Rescale transfer function
  bmagPWF.RescaleTransferFunction(0.015, 0.025)
  
  
  
  # set active source
  SetActiveSource(calculator1_1)
  # set scalar coloring
  ColorBy(calculator1_1Display, ('CELLS', 'Jmag'))
  # rescale color and/or opacity maps used to include current data range
  calculator1_1Display.RescaleTransferFunctionToDataRange(True, False)
  # show color bar/color legend
  calculator1_1Display.SetScalarBarVisibility(renderView1, True)
  # get color transfer function/color map for 'Jmag'
  jmagLUT = GetColorTransferFunction('Jmag')
  # change representation type
  calculator1_1Display.SetRepresentationType('Volume')
  # Properties modified on calculator1_1Display
  calculator1_1Display.ShowIsosurfaces = 1
  calculator1_1Display.IsosurfaceValues = [0.25, 0.31, 0.5]
  # hide color bar/color legend
  calculator1_1Display.SetScalarBarVisibility(renderView1, False)
  # Apply a preset using its name. Note this may not work as expected when presets have duplicate names.
  jmagLUT.ApplyPreset('Black, Blue and White', True)
  # get color transfer function/color map for 'Jmag'
  jmagLUT.RescaleTransferFunction(0.2, 0.6) #(0.001, 0.8)
  # get opacity transfer function/opacity map for 'Jmag'
  jmagPWF = GetOpacityTransferFunction('Jmag')
  # Rescale transfer function
  jmagPWF.RescaleTransferFunction(0.2, 0.6) # (0.001, 0.8)
  
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

  name1 = 'Bzero_J3_mzmy_'+str(a)
  name2 = 'Bzero_J3_mzmx_'+str(a)
  name3 = 'Bzero_J3_mxmy_'+str(a)

  # current camera placement for renderView1
  renderView1.CameraPosition = [-129.59156150296562, 12.00055973418057, 62.50025046057999]
  renderView1.CameraFocalPoint = [11.998267595332418, 12.00055973418057, 62.50025046057999]
  renderView1.CameraViewUp = [0.0, -1.0, 2.220446049250313e-16]
  renderView1.CameraParallelScale = 64.92088015465893
  
  # save screenshot
  #SaveScreenshot('/disk/plasmaz/jaa/Critical_Balance/CB8waves_104_1/Video_maglines/Bzero_J3_mzmy_20.png', renderView1, ImageResolution=[1550, 803]) #HERE SAVE THE PIC
  SaveScreenshot('/disk/plasmaz/jaa/Critical_Balance/Temperature_test/This_are_the_runs_2/one/'+ name1 +'.png', renderView1, ImageResolution=[1550, 803]) #  
  # reset view to fit data
  renderView1.ResetCamera()
  
  # current camera placement for renderView1
  renderView1.CameraPosition = [11.998267595332418, 153.5903888324786, 62.50025046057999]
  renderView1.CameraFocalPoint = [11.998267595332418, 12.00055973418057, 62.50025046057999]
  renderView1.CameraViewUp = [-1.0, 0.0, 2.220446049250313e-16]
  renderView1.CameraParallelScale = 64.92088015465893
  
  # save screenshot
  #SaveScreenshot('/disk/plasmaz/jaa/Critical_Balance/CB8waves_104_1/Video_maglines/Bzero_J3_mzmx_20.png', renderView1, ImageResolution=[1550, 803]) #HERE SAVE THE PIC
  SaveScreenshot('/disk/plasmaz/jaa/Critical_Balance/Temperature_test/This_are_the_runs_2/one/'+ name2 +'.png', renderView1, ImageResolution=[1550, 803]) 
  # reset view to fit data
  renderView1.ResetCamera()
  
  # current camera placement for renderView1
  renderView1.CameraPosition = [11.998267595332418, 12.00055973418057, 179.5166381451238]
  renderView1.CameraFocalPoint = [11.998267595332418, 12.00055973418057, 62.50025046057999]
  renderView1.CameraViewUp = [4.440892098500626e-16, -1.0, 0.0]
  renderView1.CameraParallelScale = 64.92088015465893
  
  # save screenshot
  #SaveScreenshot('/disk/plasmaz/jaa/Critical_Balance/CB8waves_104_1/Video_maglines/Bzero_J3_mxmy_20.png', renderView1, ImageResolution=[1550, 803]) #HERE SAVE THE PIC
  SaveScreenshot('/disk/plasmaz/jaa/Critical_Balance/Temperature_test/This_are_the_runs_2/one/'+ name3 +'.png', renderView1, ImageResolution=[1550, 803])
  #### saving camera placements for all active views
  
  # current camera placement for renderView1
  renderView1.CameraPosition = [11.998267595332418, 12.00055973418057, 179.5166381451238]
  renderView1.CameraFocalPoint = [11.998267595332418, 12.00055973418057, 62.50025046057999]
  renderView1.CameraViewUp = [4.440892098500626e-16, -1.0, 0.0]
  renderView1.CameraParallelScale = 64.92088015465893
  
#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).