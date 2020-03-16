#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  7 11:05:51 2019

@author: jaa
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 25 11:06:02 2019

@author: jaa
"""

import numpy as np
import matplotlib.pyplot as plt
import h5py
import glob, os




from mpl_toolkits import mplot3d
from matplotlib import cm
from mpl_toolkits.mplot3d import axes3d


from scipy import signal
from scipy.stats import linregress
from matplotlib.ticker import NullFormatter  # useful for `logit` scale


#cd /run/media/jaa/C2BCB9BCBCB9AB75/DIRAC_RUNs/SCALING/weak_2_scaling_256

#Call and read the file
###############################################################################


#Load the file
#for i in range(1,10):
#path = '/disk/plasma2/jaa/CB8WAVES/CB8waves_04'
os.chdir('/disk/plasmaz/jaa/Critical_Balance/CB8waves_104_1/Bxy_xytime/')
path = '/disk/plasmaz/jaa/Critical_Balance/CB8waves_104_1'

file=sorted(glob.glob(os.path.join(path, '*.h5')))
a=-1
for i in glob.glob( os.path.join(path, '*.h5')):
    a=a+1;   
#    print(file[a])   
    filename1= file[5]#'/disk/plasma2/jaa/CB8WAVES/CB8waves_04/pfd.002000_p000000.h5' 
    #Read the file
    print(filename1)
    h5_1st = h5py.File(filename1, 'r')
    #Get the fields
    ###############################################################################
    e_field=h5_1st[list(h5_1st.keys())[16]] # Like this the outcome is a group
    h_field=h5_1st[list(h5_1st.keys())[17]] # ['hx', 'hy', 'hz']
    j_field=h5_1st[list(h5_1st.keys())[18]] # ['jx', 'jy', 'jz']
    v_field=h5_1st[list(h5_1st.keys())[24]] # ['vx_e', 'vx_i', 'vy_e', 'vy_i', 'vz_e', 'vz_i']
    #n_density=h5_1st[list(h5_1st.keys())[22]] #['n_e', 'n_i']
    #T_tensor=h5_1st[list(h5_1st.keys())[0]] #['Txx_e', 'Txx_i', 'Txy_e', 'Txy_i', 'Txz_e', 'Txz_i', 'Tyy_e', 'Tyy_i', 'Tyz_e', 'Tyz_i', 'Tzz_e', 'Tzz_i']

    # Electric field components
    #-----------------------------------------------------
    
    #ex=e_field['ex'] ; ex=ex['p0'] ; ex=ex['3d'] ; ex=np.transpose(ex, (1, 2, 0))
    #ey=e_field['ey'] ; ey=ey['p0'] ; ey=ey['3d'] ; ey=np.transpose(ey, (1, 2, 0))
    #ez=e_field['ez'] ; ez=ez['p0'] ; ez=ez['3d'] ; ez=np.transpose(ez, (1, 2, 0))

    # Magnetic field components
    #-----------------------------------------------------
    hx=h_field['hx'] ; hx=hx['p0'] ; hx=hx['3d'] ; hx=np.transpose(hx, (1, 2, 0))
    hy=h_field['hy'] ; hy=hy['p0'] ; hy=hy['3d'] ; hy=np.transpose(hy, (1, 2, 0))
    hz=h_field['hz'] ; hz=hz['p0'] ; hz=hz['3d'] ; hz=np.transpose(hz, (1, 2, 0))

    # Current field components
    #-----------------------------------------------------
    jx=j_field['jx'] ; jx=jx['p0'] ; jx=jx['3d'] ; jx=np.transpose(jx, (1, 2, 0))
    jy=j_field['jy'] ; jy=jy['p0'] ; jy=jy['3d'] ; jy=np.transpose(jy, (1, 2, 0))
    jz=j_field['jz'] ; jz=jz['p0'] ; jz=jz['3d'] ; jz=np.transpose(jz, (1, 2, 0))

    # Velocity Ions components
    #-----------------------------------------------------
    vx_i=v_field['vx_i'] ; vx_i=vx_i['p0'] ; vx_i=vx_i['3d'] ; vx_i=np.transpose(vx_i, (1, 2, 0))
    vy_i=v_field['vy_i'] ; vy_i=vy_i['p0'] ; vy_i=vy_i['3d'] ; vy_i=np.transpose(vy_i, (1, 2, 0))
    vz_i=v_field['vz_i'] ; vz_i=vz_i['p0'] ; vz_i=vz_i['3d'] ; vz_i=np.transpose(vz_i, (1, 2, 0))

    # Velocity electrons components
    #-----------------------------------------------------
    vx_e=v_field['vx_e'] ; vx_e=vx_e['p0'] ; vx_e=vx_e['3d'] ; vx_e=np.transpose(vx_e, (1, 2, 0))
    vy_e=v_field['vy_e'] ; vy_e=vy_e['p0'] ; vy_e=vy_e['3d'] ; vy_e=np.transpose(vy_e, (1, 2, 0))
    vz_e=v_field['vz_e'] ; vz_e=vz_e['p0'] ; vz_e=vz_e['3d'] ; vz_e=np.transpose(vz_e, (1, 2, 0))

    # Density Ions and electrons
    #-----------------------------------------------------
    #n_i=n_density['n_i'] ; n_i=n_i['p0'] ; n_i=n_i['3d'] ; n_i=np.transpose(n_i, (1, 2, 0))
    #n_e=n_density['n_e'] ; n_e=n_e['p0'] ; n_e=n_e['3d'] ; n_e=np.transpose(n_e, (1, 2, 0))

    """
    # Temperature electron tensor components
    #-----------------------------------------------------
    Txx_e = T_tensor['Txx_e'] ; Txx_e=Txx_e['p0'] ; Txx_e=Txx_e['3d'] ; Txx_e=np.transpose(Txx_e, (1, 2, 0))
    Txy_e = T_tensor['Txy_e'] ; Txy_e=Txy_e['p0'] ; Txy_e=Txy_e['3d'] ; Txy_e=np.transpose(Txy_e, (1, 2, 0))
    Txz_e = T_tensor['Txz_e'] ; Txz_e=Txz_e['p0'] ; Txz_e=Txz_e['3d'] ; Txz_e=np.transpose(Txz_e, (1, 2, 0))
    Tyy_e = T_tensor['Tyy_e'] ; Tyy_e=Tyy_e['p0'] ; Tyy_e=Tyy_e['3d'] ; Tyy_e=np.transpose(Tyy_e, (1, 2, 0))
    Tyz_e = T_tensor['Tyz_e'] ; Tyz_e=Tyz_e['p0'] ; Tyz_e=Tyz_e['3d'] ; Tyz_e=np.transpose(Tyz_e, (1, 2, 0))
    Tzz_e = T_tensor['Tzz_e'] ; Tzz_e=Tzz_e['p0'] ; Tzz_e=Tzz_e['3d'] ; Tzz_e=np.transpose(Tzz_e, (1, 2, 0))
    # Temperature Ions tensor components
    #-----------------------------------------------------
    Txx_i = T_tensor['Txx_i'] ; Txx_i=Txx_i['p0'] ; Txx_i=Txx_i['3d'] ; Txx_i=np.transpose(Txx_i, (1, 2, 0))
    Txy_i = T_tensor['Txy_i'] ; Txy_i=Txy_i['p0'] ; Txy_i=Txy_i['3d'] ; Txy_i=np.transpose(Txy_i, (1, 2, 0))
    Txz_i = T_tensor['Txz_i'] ; Txz_i=Txz_i['p0'] ; Txz_i=Txz_i['3d'] ; Txz_i=np.transpose(Txz_i, (1, 2, 0))
    Tyy_i = T_tensor['Tyy_i'] ; Tyy_i=Tyy_i['p0'] ; Tyy_i=Tyy_i['3d'] ; Tyy_i=np.transpose(Tyy_i, (1, 2, 0))
    Tyz_i = T_tensor['Tyz_i'] ; Tyz_i=Tyz_i['p0'] ; Tyz_i=Tyz_i['3d'] ; Tyz_i=np.transpose(Tyz_i, (1, 2, 0))
    Tzz_i = T_tensor['Tzz_i'] ; Tzz_i=Tzz_i['p0'] ; Tzz_i=Tzz_i['3d'] ; Tzz_i=np.transpose(Tzz_i, (1, 2, 0))
    ###############################################################################
    """
    
  
    #Make operations with the loaded data
    ###############################################################################
    
    #plt.pcolor(np.sum(ex,axis=0)) #integration over the first dimension
    # ex.shape returns the dimensions of the matrix
    Bxy = np.sqrt(hx[:,:,:]**2 + hy[:,:,:]**2)
    a=5
    nn= int(a)
    nstr=str(nn)
    fig1, ax0 = plt.subplots()
    im = plt.pcolor(Bxy[:,0,:]) #plot in the 65 (middle)
    fig1.colorbar(im, ax=ax0)
    ax0.set_title('Bxy')
    fig2, ax0 = plt.subplots()
    im = plt.pcolor(Bxy[0,:,:]) #plot in the 65 (middle)
    fig2.colorbar(im, ax=ax0)
    ax0.set_title('Bxy')
    fig1.savefig('Bxy_xz'+nstr+'.png', bbox_inches='tight')
    fig2.savefig('Bxy_yz'+nstr+'.png', bbox_inches='tight')
#Make the plot
#####################################################################



    
    #3D Contour plot
    ###############################################################################
    
    
    
    fig= plt.figure()
    ax = plt.axes(projection='3d')
    ax.contour3D(Bxy[5,:,:], Bxy[5,:,:], Bxy[5,:,:], 50, cmap='viridis')
    ax.set_xlabel('hx')
    ax.set_ylabel('hy')
    ax.set_zlabel('hz')
    ax.view_init(60, 35)
    fig
    #------------------------------------------------------------------------
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.plot_wireframe(hx[5,:,:], hy[5,:,:], hz[5,:,:], color='black')
    ax.set_title('wireframe');
    #------------------------------------------------------------------------
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.plot_surface(hx[5,:,:], hy[5,:,:], hz[5,:,:], rstride=1, cstride=1,
                    cmap='viridis', edgecolor='none')
    ax.set_title('surface');
    ###############################################################################
    
    fig, ax0 = plt.subplots()
    im = plt.pcolor(ex[65,:,:]) #plot in the 65 (middle)
    fig.colorbar(im, ax=ax0)
    ax0.set_title('pcolormesh with levels')
    


X=np.arange(400)
Y=np.arange(400)
Z=np.arange(2016)

mesh = np.meshgrid(X,Y, Z)

xyv, yxv= np.meshgrid(X, Y)
xzv, zxv = np.meshgrid(Z, X)
yzv, yzv = np.meshgrid(Z, Y)


fig = plt.figure()
ax = fig.gca(projection='3d')
X, Y, Z = axes3d.get_test_data(0.05)
#ax.plot_surface(X, Y, Z, rstride=8, cstride=8, alpha=0.3)
cset = ax.contour(xzv, xzv, xzv, zdir='z', offset=-100, cmap=cm.coolwarm)
#cset = ax.contour(X, Y, Z, zdir='x', offset=-40, cmap=cm.coolwarm)
#cset = ax.contour(X, Y, Z, zdir='y', offset=40, cmap=cm.coolwarm)



plt.show()
ax.set_xlabel('X')
ax.set_xlim(-40, 40)
ax.set_ylabel('Y')
ax.set_ylim(-40, 40)
ax.set_zlabel('Z')
ax.set_zlim(-100, 100)

#cd /disk/plasma2/jaa/CB8WAVES/CB8waves_04/



fig = go.Figure(data=[
    go.Mesh3d(
        x=X,
        y=Y,
        z=Z,
        colorbar_title='z',
        colorscale=[[0, 'gold'],
                    [0.5, 'mediumturquoise'],
                    [1, 'magenta']],
        # Intensity of each vertex, which will be interpolated and color-coded
        intensity=[0, 0.33, 0.66, 1],
        # i, j and k give the vertices of triangles
        # here we represent the 4 triangles of the tetrahedron surface
        i=[0, 0, 0, 1],
        j=[1, 2, 3, 2],
        k=[2, 3, 1, 3],
        name='y',
        showscale=True
    )
])

fig.show()


