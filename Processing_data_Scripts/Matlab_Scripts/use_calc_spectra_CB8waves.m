%#-----------------------------------------------------------------------
%C= 10^{-3}; beta =1
%#----------------------------------------------------------------------
%diary myDiaryFile
%cd '/disk/plasma2/jaa/CB8WAVES/CB8waves_04'; %Set the directory of the files 
%path = '/disk/plasma2/jaa/CB8WAVES/CB8waves_04';

cd '/disk/plasmaz/jaa/Critical_Balance/CB8waves_104_1/Raw_outputs'; %This is important because the xdmf files are in that directory
path = '/disk/plasmaz/jaa/Critical_Balance/CB8waves_104_1/Raw_outputs';

S = dir(fullfile(path,'*.h5'));
N=numel(S); % number of files to use
H=zeros(N);

for i=2:N
    disp(strcat('Computing step ...',S(i).name))
    
    fileID =  S(i).name; %change the 1 per i
    %h5disp(fileID);
    info = h5info(fileID);
    %T_1st = info.Groups(1).Name;
    %E_1st = info.Groups(17).Name;
    B_1st = info.Groups(18).Name;
    J_1st = info.Groups(19).Name;
    %V_1st = info.Groups(25).Name;
    %n_1st = info.Groups(23).Name;
    %diary off
    %s2='/hx/p0/3d';
    %sx = strcat(B_1st,s2);
    Bx=h5read(fileID,strcat(B_1st,'/hx/p0/3d'));
    By=h5read(fileID,strcat(B_1st,'/hy/p0/3d'));
    Bz=h5read(fileID,strcat(B_1st,'/hz/p0/3d'));
    %Ex=h5read(fileID,strcat(E_1st,'/ex/p0/3d'));
    %Ey=h5read(fileID,strcat(E_1st,'/ey/p0/3d'));
    %Ez=h5read(fileID,strcat(E_1st,'/ez/p0/3d'));
    %vix=h5read(fileID,strcat(V_1st,'/vx_i/p0/3d'));
    %viy=h5read(fileID,strcat(V_1st,'/vy_i/p0/3d'));
    %viz=h5read(fileID,strcat(V_1st,'/vz_i/p0/3d'));
    %vex=h5read(fileID,strcat(V_1st,'/vx_e/p0/3d'));
    %vey=h5read(fileID,strcat(V_1st,'/vy_e/p0/3d'));
    %vez=h5read(fileID,strcat(V_1st,'/vz_e/p0/3d'));
    %ni = h5read(fileID,strcat(n_1st,'/n_i/p0/3d'));
    %ne = h5read(fileID,strcat(n_1st,'/n_e/p0/3d'));
    
    Jx=h5read(fileID,strcat(J_1st,'/jx/p0/3d'));
    Jy=h5read(fileID,strcat(J_1st,'/jy/p0/3d'));
    Jz=h5read(fileID,strcat(J_1st,'/jz/p0/3d'));
    
    %------------------------------------------------------------------
    %resx=18/296; %0.06081
    %resy=18/296;
    %resz=63/1024; %0.06152
    
    res=0.06;
    threshold=52; %100; % Check if by changing this the extension in the kper will be reduced
    
    %for 103
    %res=0.07;
    
    cd '/disk/plasma2/jaa/CB8WAVES/Processing_data_Scripts/Matlab_Scripts/'
    spectraB=CalcSpectra_CB(Bx,By,Bz,res,threshold);
    %spectraE=CalcSpectra_CB(Ex,Ey,Ez,res,threshold);
    %spectravi=CalcSpectra_CB(vix,viy,viz,res,threshold);
    %spectrave=CalcSpectra_CB(vex,vey,vez,res,threshold);
    %spectrani=CalcSpectra_CB_n(ni,res,threshold); %here this is a different file
    %spectrane=CalcSpectra_CB_n(ne,res,threshold);
    spectraJ=CalcSpectra_CB(Jx,Jy,Jz,res,threshold);
    %cd '/disk/plasma2/jaa/CB8WAVES/CB8waves_03' %I must bear this in mind
    %cd '/disk/plasmaz/jaa/Critical_Balance/CB8waves_104_1/';
    
    
    
    % Making plots
    %----------------------------------------------------------------------
    
    
    f1=figure(1);
    [C,h]=contourf(spectraB.kpar,spectraB.kper,log10(spectraB.P2D'),'LevelStep',0.5);
    %[C]=contourf(spectraB.kper,spectraB.kpar,log10(spectraE.P2D),'LevelStep',0.5);
    set(h,'LineColor','none')
    lim = caxis; 
    caxis([0 13]);
    colormap(jet)
    colorbar
    set(gca,'XScale','log','YScale','log')
    xlabel('$k_\| d_i$','Interpreter','latex')
    ylabel('$k_\perp d_i$','Interpreter','latex')
    title('B','Interpreter','latex')
    x=log10([0.1 1]);
    x3=log10([1 70]);
    y=(3/2)*x;
    y3=(3)*x3;
    hold on
    loglog(10.^x,10.^y,'--k')
    loglog(10.^x3,10.^y3,'--b')
    ylim([0.28 52])
    txt = '\bf{3/2 \rightarrow}';
    txt2 = '\bf{\leftarrow 3}';
    text(0.2,0.5,txt)
    text(5,28,txt2)
    hold off
    
    %{
    f12=figure(12);
    [C,h]=contourf(spectraB.kper,spectraB.kpar,log10(spectraE.P2D),'LevelStep',0.5);
    set(h,'LineColor','none')
    colormap(hot)
    colorbar
    set(gca,'XScale','log','YScale','log')
    ylabel('$k_\| d_i$','Interpreter','latex')
    xlabel('$k_\perp d_i$','Interpreter','latex')
    title('B','Interpreter','latex')
    x=log10([0.26 1]);
    x3=log10([1 70]);
    y=(2/3)*(x-0.314)+(1/5);
    y3=(1/3)*(x3-1)+(1/3);
    %y1=get(gca,'ylim'); %x1=10; %x2=get(gca,'xlim'); y2=10;
    %grid on
    hold on
    loglog(10.^(x),10.^y,'--k')
    loglog(10.^(x3),10.^y3,'--b')
    %loglog([x1 x1],y1,'--k') %loglog(x2,y2*ones(size(x)),'--k')
    %pbaspect([1 1 1])
    %%xlim([0.1 50])
    ylim([0.28 50])
    %legend('2/3','1/2')
    txt = '\bf{\leftarrow 2/3}';
    txt2 = '\bf{\downarrow 1/3}';
    text(0.8,0.7,txt)
    text(20,4,txt2)
    hold off
       
    
    f2=figure(2);
    [C,h]=contourf(spectraE.kpar,spectraE.kper,log10(spectraE.P2D'),'LevelStep',0.5);
    set(h,'LineColor','none')
    colormap(hot)
    colorbar
    set(gca,'XScale','log','YScale','log')
    xlabel('$k_\| d_i$','Interpreter','latex')
    ylabel('$k_\perp d_i$','Interpreter','latex')
    title('E','Interpreter','latex')
    x=log10([0.1 1]);
    x3=log10([1 70]);
    y=(3/2)*x;
    y3=(3)*x3;
    hold on
    loglog(10.^x,10.^y,'--k')
    loglog(10.^x3,10.^y3,'--b')
    ylim([0.28 50])
    txt = '\bf{3/2 \rightarrow}';
    txt2 = '\bf{\leftarrow 3}';
    text(0.2,0.5,txt)
    text(5,28,txt2)
    hold off
    %{
    f31=figure(31);
    [C,h]=contourf(spectravi.kpar,spectravi.kper,log10(spectravi.P2D'),'LevelStep',0.5);
    set(h,'LineColor','none')
    colorbar
    colormap(hot)
    set(gca,'XScale','log','YScale','log')
    xlabel('$k_\| d_i$','Interpreter','latex')
    ylabel('$k_\perp d_i$','Interpreter','latex')
    title('$v_i$','Interpreter','latex')
    ylim([1 70])
    pbaspect([1 1 1])
    %}
    f3=figure(3);
    [C,h]=contourf(spectravi.kpar,spectravi.kper,log10(spectravi.P2D'),'LevelStep',0.5);
    set(h,'LineColor','none')
    colormap(hot)
    colorbar
    set(gca,'XScale','log','YScale','log')
    xlabel('$k_\| d_i$','Interpreter','latex')
    ylabel('$k_\perp d_i$','Interpreter','latex')
    title('$v_i$','Interpreter','latex')
    x=log10([0.1 1]);
    x3=log10([1 70]);
    y=(3/2)*x;
    y3=(3)*x3;
    hold on
    loglog(10.^x,10.^y,'--k')
    loglog(10.^x3,10.^y3,'--b')
    ylim([0.28 50])
    txt = '\bf{3/2 \rightarrow}';
    txt2 = '\bf{\leftarrow 3}';
    text(0.2,0.5,txt)
    text(5,28,txt2)
    %ylim([1 70])
    %pbaspect([1 1 1])
    hold off
    
    
    f4=figure(4);
    [C,h]=contourf(spectrave.kpar,spectrave.kper,log10(spectrave.P2D'),'LevelStep',0.5);
    set(h,'LineColor','none')
    colormap(hot)
    colorbar
    set(gca,'XScale','log','YScale','log')
    xlabel('$k_\| d_i$','Interpreter','latex')
    ylabel('$k_\perp d_i$','Interpreter','latex')
    title('$v_e$','Interpreter','latex')
    x=log10([0.1 1]);
    x3=log10([1 70]);
    y=(3/2)*x;
    y3=(3)*x3;
    hold on
    loglog(10.^x,10.^y,'--k')
    loglog(10.^x3,10.^y3,'--b')
    ylim([0.28 50])
    txt = '\bf{3/2 \rightarrow}';
    txt2 = '\bf{\leftarrow 3}';
    text(0.2,0.5,txt)
    text(5,28,txt2)
    hold off
    
    f41=figure(41);
    [C,h]=contourf(spectrani.kpar,spectrani.kper,log10(spectrani.P2D'),'LevelStep',0.5);
    set(h,'LineColor','none')
    colormap(hot)
    colorbar
    set(gca,'XScale','log','YScale','log')
    xlabel('$k_\| d_i$','Interpreter','latex')
    ylabel('$k_\perp d_i$','Interpreter','latex')
    title('$n_i$','Interpreter','latex')
    x=log10([0.1 1]);
    x3=log10([1 70]);
    y=(3/2)*x;
    y3=(3)*x3;
    hold on
    loglog(10.^x,10.^y,'--k')
    loglog(10.^x3,10.^y3,'--b')
    ylim([0.28 50])
    txt = '\bf{3/2 \rightarrow}';
    txt2 = '\bf{\leftarrow 3}';
    text(0.2,0.5,txt)
    text(5,28,txt2)
    hold off
    %}
    f5=figure(5);
    %loglog(spectraB.kper,spectraB.P1Dper/rms(spectraB.P1Dper), 'k', spectraE.kper,spectraE.P1Dper/rms(spectraE.P1Dper), 'b',spectravi.kper,spectravi.P1Dper/rms(spectravi.P1Dper), 'r', spectrave.kper,spectrave.P1Dper/rms(spectrave.P1Dper), 'g', spectrani.kper,spectrani.P1Dper/rms(spectrani.P1Dper), 'c',spectrane.kper,spectrane.P1Dper/rms(spectrane.P1Dper), 'm')
    %loglog(spectraB.kper,spectraB.P1Dper/rms(spectraB.P1Dper), 'k',spectravi.kper,spectravi.P1Dper/rms(spectravi.P1Dper), 'r')
    loglog(spectraB.kper,spectraB.P1Dper/rms(spectraB.P1Dper), 'k',spectraJ.kper,spectraJ.P1Dper/rms(spectraJ.P1Dper), 'b')
    %x2=log10([1 16]); y2=x2.*(-3.2)+14.3; x3=log10([16 50]); y3=x3.*(-0.3)+10.8;
    x2=log10([1 10]); y2=x2.*(-2.8); x3=log10([0.099 1]); y3=x3.*(-1.7);
    x33=1;    x44=10;    x55=14.2857;
    y1=get(gca,'ylim');
    hold on
    loglog(10.^x2,10.^y2,'--k')
    loglog(10.^x3,10.^y3,'--b')
    loglog([x33 x33],y1,'--')
    loglog([x44 x44],y1,'--')
    loglog([x55 x55],y1,'--')
    % Text in the plot
    txt3 = '$$ k_{\perp}d_e=1$$';
    txt4 = '$$ k_{\perp}d_i=1$$';
    txt5 = '$$ k_{\perp}^{-2.8}$$';
    txt6 = '$$ k_{\perp}^{-1.7}$$';
    txt8 = '$$ k_{\perp}\lambda_D=1$$';
    text(5 , 30,txt3,'Interpreter','latex'); text(1.1 , 30,txt4,'Interpreter','latex'); text(1.1 , 0.1,txt5,'Interpreter','latex')
    text(0.35 , 16,txt6,'Interpreter','latex'); text(16 , 0.3,txt8,'Interpreter','latex')
    %legend({'$\tilde{B}$','$\tilde{E}$','$\tilde{v}_i$','$\tilde{v}_e$','$\tilde{n}_i$', '$\tilde{n}_e$'},'Interpreter','latex')
    %legend({'$\tilde{B}$','$\tilde{v}_i$'},'Interpreter','latex')
    legend({'$\tilde{B}$','$\tilde{J}$'},'Interpreter','latex')
    %legend({'$B$','$E$','$v_i$','$v_e$','$n_i$'},'Interpreter','latex')
    xlabel('$k_\perp d_i$','Interpreter','latex')
    ylabel('$P1D_\perp$','Interpreter','latex')
    xlim([0.3 52.36]); ylim([10^(-5) 100])
    title('$P1D_\perp$','Interpreter','latex')
    hold off
    
    
    f6=figure(6);
    %loglog(spectraB.kpar,spectraB.P1Dpar/rms(spectraB.P1Dpar), 'k', spectraE.kpar,spectraE.P1Dpar/rms(spectraE.P1Dpar), 'b',spectravi.kpar,spectravi.P1Dpar/rms(spectravi.P1Dpar), 'r', spectrave.kpar,spectrave.P1Dpar/rms(spectrave.P1Dpar), 'g', spectrani.kpar,spectrani.P1Dpar/rms(spectrani.P1Dpar), 'c',spectrane.kpar,spectrane.P1Dpar/rms(spectrane.P1Dpar), 'm')
    %loglog(spectraB.kpar,spectraB.P1Dpar/rms(spectraB.P1Dpar), 'k',spectravi.kpar,spectravi.P1Dpar/rms(spectravi.P1Dpar), 'r')
    loglog(spectraB.kpar,spectraB.P1Dpar/rms(spectraB.P1Dpar), 'k',spectraJ.kpar,spectraJ.P1Dpar/rms(spectraJ.P1Dpar), 'b')
    x=log10([0.099 1]);
    x2=log10([1 10]);
    x3=log10([5 30]);
    %y=x.*(-2.3)+12.2;
    y=x.*(-2);%+12.2;
    y2=x2.*(-3.5);%+12.2;
    y3=x3.*(-0.3);%+10.6;
    x33=1;    x44=10;   x55=14.2857;
    xlim([0.3 52.36]); ylim([10^(-5) 100]);
    y1=get(gca,'ylim');
    hold on
    loglog(10.^x,10.^y,'--b')
    loglog(10.^x2,10.^y2,'--k')
    %loglog(10.^x3,10.^y3,'--')
    loglog([x33 x33],y1,'--')
    loglog([x44 x44],y1,'--')
    loglog([x55 x55],y1,'--')
    % Text in the plot
    xlabel('$k_\| d_i$','Interpreter','latex')
    ylabel('$P1D_\|$','Interpreter','latex')
    txt3 = '$$ k_{||}d_e=1$$'; txt4 = '$$ k_{||}d_i=1$$';    txt5 = '$$ k_{||}^{-3.5}$$';
    txt6 = '$$ k_{||}^{-0.3}$$';    txt7 = '$$ k_{||}^{-2}$$';    txt8 = '$$ k_{||}\lambda_D=1$$';
    text(5 , 30,txt3,'Interpreter','latex'); text(1.1 , 30,txt4,'Interpreter','latex')
    text(2 , 0.4,txt5,'Interpreter','latex'); %text(18 , 1,txt6,'Interpreter','latex')
    text(0.5 , 7,txt7,'Interpreter','latex'); text(16 , 0.3,txt8,'Interpreter','latex')
    title('$P1D_\|$','Interpreter','latex')
    %legend({'$B$','$E$','$v_i$','$v_e$','$n_i$'},'Interpreter','latex') 
    %legend({'$\tilde{B}$','$\tilde{v}_i$'},'Interpreter','latex')
    legend({'$\tilde{B}$','$\tilde{J}$'},'Interpreter','latex')
    hold off
    
    %To plot rms values and also make large the text
    %{ 
    f111=figure(111);
    plot(Rmsvalues.t_pi,Rmsvalues.B_rms);
    hold on;
    plot(Rmsvalues.t_pi,Rmsvalues.E_rms);
    plot(Rmsvalues.t_pi,Rmsvalues.J_rms);
    xlabel('$t \omega_{pi}$','Interpreter','latex')
    %ylabel('$P1D_\|$','Interpreter','latex')
    legend({'$B_{rms}$','$E_{rms}$','$J_{rms}$'},'Interpreter','latex')
    set(gca,'FontSize',15)
    hold off;    
    saveas(f111,'RMS_values.png');
    %}
    
    
    %----------------------------------------------------------------------
    cd /disk/plasmaz/jaa/Critical_Balance/CB8waves_104_1/Statistics/2D_1D_PSD
    % Save the plots
    %-------------------------------------------------------------------------
    %saveas(f1,strcat(S(i).name,'2Dspectrum_B_CB04.png'));
    %saveas(f12,strcat(S(i).name,'2Dspectrum_B_trans_CB04.png'));
    %saveas(f2,strcat(S(i).name,'2Dspectrum_E_CB04.png'));
    %saveas(f3,strcat(S(i).name,'2Dspectrum_vi_CB04.png'));
    %saveas(f4,strcat(S(i).name,'2Dspectrum_ve_CB04.png'));
    %saveas(f41,strcat(S(i).name,'2Dspectrum_ni_CB04.png'));
    saveas(f5,strcat(S(i).name,'1_perDspectrum_CB04_BJ.png'));
    saveas(f6,strcat(S(i).name,'1_parDspectrum_CB04_BJ.png'));
    cd '/disk/plasmaz/jaa/Critical_Balance/CB8waves_104_1/Raw_outputs'; 
    
end
