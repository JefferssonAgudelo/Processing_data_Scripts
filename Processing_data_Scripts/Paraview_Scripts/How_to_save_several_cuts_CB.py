# trace generated using paraview version 5.5.0

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

pfd0 = FindSource('pfd.0*')
slice1 = FindSource('Slice1')
calculator2 = FindSource('Calculator2')
pfd0_1 = FindSource('pfd.0*')
calculator1 = FindSource('Calculator1')
slice2 = GetActiveSource()


for i in range(20):
    a=float(i+1)
    slice2.SliceType.Origin = [10.0, 10.0, a]
    renderView2 = GetActiveViewOrCreate('RenderView')
    renderView2.Update()
    renderView2.ResetCamera()
    renderView2.CameraPosition = [10.0, 10.0, 56.64101615137755]
    renderView2.CameraFocalPoint = [10.0, 10.0, a]
    renderView2.CameraParallelScale = 14.142135623730951
    SaveScreenshot('/run/media/jaa/C2BCB9BCBCB9AB75/DIRAC_RUNs/SCALING/256_1_mode/New Folder/jz_600_' + str(int(a)) + '.png', renderView2, ImageResolution=[1083, 535])
    renderView2.CameraPosition = [10.0, 10.0, 56.64101615137755]
    renderView2.CameraFocalPoint = [10.0, 10.0, a]
    renderView2.CameraParallelScale = 14.142135623730951


slice2 = GetActiveSource()
for i in range(20):
    a=float(i+1)
    slice2.SliceType.Origin = [a, 10.0, 10.0]
    renderView2 = GetActiveViewOrCreate('RenderView')
    renderView2.Update()
    renderView2.ResetCamera()
    renderView2.CameraPosition = [56.64101615137755, 10.0, 10.0]
    renderView2.CameraFocalPoint = [a, 10.0, 10.0]
    renderView2.CameraParallelScale = 14.142135623730951
    a=int(a)
    SaveScreenshot('/run/media/jaa/C2BCB9BCBCB9AB75/DIRAC_RUNs/SCALING/256_1_mode/New Folder/jz_600_x_' + str(int(a)) + '.png', renderView2, ImageResolution=[1083, 535])
    renderView2.CameraPosition = [56.64101615137755, 10.0, 10.0]
    renderView2.CameraFocalPoint = [a, 10.0, 10.0]
    renderView2.CameraParallelScale = 14.142135623730951


#### saving camera placements for all active views

# current camera placement for renderView2


#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).