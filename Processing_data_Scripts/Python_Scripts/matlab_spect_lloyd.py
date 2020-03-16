#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 15 11:13:56 2019

@author: jaa
"""

cd '/run/media/jaa/C2BCB9BCBCB9AB75/DIRAC_RUNs/SCALING/256_1_mode'

h5disp('pfd.000600_p000000.h5')

Bx=h5read('pfd.000600_p000000.h5','/h-uid-0x11b7960/hx/p0/3d');
By=h5read('pfd.000600_p000000.h5','/h-uid-0x11b7960/hy/p0/3d');
Bz=h5read('pfd.000600_p000000.h5','/h-uid-0x11b7960/hz/p0/3d');


h5disp('pfd.030000_p000000.h5')


res=0.078125;
threshold=1000;

spectra=CalcSpectra(Bx,By,Bz,res,threshold);


[C]=contourf(spectra.kpar,spectra.kper,log10(spectra.P2D'),'LevelStep',0.5);
colorbar
set(gca,'XScale','log','YScale','log')
xlabel('$k_\| d_i$','Interpreter','latex')
ylabel('$k_\perp d_i$','Interpreter','latex')

loglog(spectra.kper,spectra.P1Dper)
x=log10([0.1 5]);
y=x.*(-5/3)+13;
x2=log10([1 15]);
y2=x2.*(-2.8)+12;
hold on
loglog(10.^x,10.^y)
loglog(10.^x2,10.^y2)
hold off
