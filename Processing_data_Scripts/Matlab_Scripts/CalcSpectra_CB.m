function spectra=CalcSpectra_CB(varx,vary,varz,res,threshold)
%This code calculates the 2D and 1D reduced power spectrum of the magnetic
%field. Bx, By, and Bz are the 3D arrays; resx,resy, resz are the spatial resolution in
%terms of d_i ; threshold is the k_perp at which I want to filter (taking the
%isocontour of the 2D spectrum corrisponding at P(k_par=0,k_per=threshold)
%This code is based on Luca Franci's routines (see Franci et al. ApJ 2018)
%and adapted by Lloyd Woodham [Current version: 15/08/2019].
%Adapted by Jeffersson Agudelo [version: 21/10/19]
disp('Computing 2D and 1D spectra...')
[kpar,kper,P2D]=Compute2Dspectrum(res,varx,vary,varz);
[P1Dper,P1Dpar]=Compute1DSpectra(P2D,kper,threshold);
spectra.kpar=kpar; spectra.kper=kper;
spectra.P2D=P2D; spectra.P1Dpar=P1Dpar; spectra.P1Dper=P1Dper;
clear varx vary varz;
end

function [kpar,kper,P2D,Nkper]=Compute2Dspectrum(res,var1,var2,var3)
%Compute the 2D power swpectrum of variables 'var1', 'var2' & 'var3'.

dx=res; %Fraction of d_i, assumes resolution is same in x/y/z directions
dy=res;
dz=res;

N=size(var1); %all the variables has the same size
Lx=N(1)*dx; %Length of box as function of d_i
Ly=N(2)*dy; 
Lz=N(3)*dz; 
    
if nargin==2 %Calculate 3D Fourier transform of each variable (narging return the number of inputs of a function)
    [ff,kx,ky,kz]=Four3D(var1,Lx,Ly,Lz);
    ff=ff.^2;
elseif nargin==3
    [fx,kx,ky,kz]=Four3D(var1,Lx,Ly,Lz);
    [fy,kx,ky,kz]=Four3D(var2,Lx,Ly,Lz);
    ff=fx.^2+fy.^2;
elseif nargin==4
    [fx,kx,ky,kz]=Four3D(var1,Lx,Ly,Lz); %3D transform of each component of B
    [fy,kx,ky,kz]=Four3D(var2,Lx,Ly,Lz);
    [fz,kx,ky,kz]=Four3D(var3,Lx,Ly,Lz);
    ff=fx.^2+fy.^2+fz.^2; %Trace
end

P3D=ff./(N(1)*N(2)*N(3)); %Normalise by 1/Nx*Ny*Nz  %Second normalization? jeff

for j=1:size(ff,2)
    for i=1:size(ff,1) 
        kp(i,j)=sqrt(ky(j).^2+kx(i).^2); %Calculate kperp array with negative wavenumbers
    end
end

for k=1:size(ff,3)
    kp2(k)=abs(kz(k)); %Calculate kpar vector with negative wavenumbers
end

%Calculate kpar and kperp for 2D spectrum:
dkpar=kz(2)-kz(1);
nkpar=fix(max(max(kz)/dkpar))+1;
kpar=(0:nkpar-1).*dkpar;

dkper=kx(2)-kx(1);
nkper=fix(max(max(kp)/dkper))+1;
kper=(0:nkper-1).*dkper;

%Use linear interpolation to intergrate in concentric circles of k_perp:

%Initialise arrays:
Ptemp=zeros(size(ff,3),nkper); %old kper vs. new kper
P2D=zeros(nkpar,nkper); %kpar vs. kper

Nkper=zeros(size(ff,3),nkper);

i1=fix(kp/dkper); %Index linking kp with each kper value
i2=i1+1;
i3=fix(kp2/dkpar); %Index linking kp2 with each kpar value

k2=kp/dkper-i1;
k1=1-k2;
    
for k= 1:size(P3D,3) %For each value of kpar, integrate over rings of kper
    P1=k1.*P3D(:,:,k); %For interpolation
    P2=k2.*P3D(:,:,k);
    for j=1:size(P3D,2) %Cycle through every grid point in x-y plane
        for i=1:size(P3D,1)
            if(i2(i,j) < nkper)
                %Use linear interpolation to sum up power to surrounding
                %rings of kper:
                Ptemp(k,i2(i,j)+1)=Ptemp(k,i2(i,j)+1)+P2(i,j); %Interpolate power to next ring
                Ptemp(k,i1(i,j)+1)=Ptemp(k,i1(i,j)+1)+P1(i,j); %Interpolate power to previous ring
                
                %Update number of points interpolated for each ring:
                Nkper(i2(i,j)+1)=Nkper(i2(i,j)+1)+k2(i,j); 
                Nkper(i1(i,j)+1)=Nkper(i1(i,j)+1)+k1(i,j);                
             end
        end
    end
    if i3(k)<nkpar
        P2D(i3(k)+1,:)=P2D(i3(k)+1,:)+Ptemp(k,:); %Fill array with kper values for each kpar
    end
end

P2D=P2D./kper(2); %Normalise // P2D=P2D./(kper(2)*kpar(2));
end

function [ft,kx,ky,kz]=Four3D(var,Lx,Ly,Lz)
    %Calculates 3D spatial Fourier transfrom of variable 'var', shifted so
    %that the centre is the zero k. xx/yy/zz are Lx, Ly and Lz.
    
    var=var-mean(mean(mean((var)))); %Detrend mean

    N=size(var); %N=[Nx,Ny,Nz]

    f=abs(fftn(var)).*sqrt(N(1)*N(2)*N(3)); %Take abs value and normalise
    %^Should it be fftn instead of ifftn?
    
    ft=fftshift(f); %Shift ft to centre wavenumber
    
    kx=((0:N(1)-1)-N(1)/2).*2*pi/Lx; %Includes negative wavenumbers
    ky=((0:N(2)-1)-N(2)/2).*2*pi/Ly;
    kz=((0:N(3)-1)-N(3)/2).*2*pi/Lz;
end

function [P1Dper,P1Dpar]=Compute1DSpectra(P2D,kper,threshold)
%Calculate parallel and perpendicular 1D sepctra from the 2D spectrum.
    [c,index]=min(abs(kper-threshold)); %Cut spectrum above a threshold k
    Pthres=P2D(2,index);    
    mask2D=P2D>Pthres;
    fieldMASKED=P2D.*mask2D;
    P1Dper=sum(fieldMASKED,1); %Sum over parallel k
    P1Dpar=sum(fieldMASKED,2); %Sum over perpendicular k
end