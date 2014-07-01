% Creates test HAWT geometry file
% Add geom creation scripts to path
addpath('./CreateGeom/')
global   ix R  CRr  bTwist HubRR Tilt eta bCone bi NBlade NElem;


% Output filename
FN=['New_V27' num2str(ix) '.geom'];

% Plot data?
PlotTurbine=0;

% Convert
dToR=pi/180;

% Create basic HAWT
Type='HAWT';
T=CreateTurbine(NBlade,NElem,0,0,R,[],[],[],Type,1,HubRR,CRr,bTwist,bi,eta,bCone,Tilt);

% Write geom file
WriteTurbineGeom(FN,T);

% Plot if desired
if PlotTurbine
    
    % Plot animated turbine rotation
    XLim=[-1,1];
    YLim=[-1,1];
    ZLim=[-1,1];
    
    % Plot controls
    PlotVec=1;
    SFVec=.5;
    Trans=.5;
    
    hf=figure(1);
    set(hf,'Position',[303   124   956   610]) 
    set(gca,'Position',[5.2743e-002  5.1245e-002  8.9979e-001  8.8141e-001])
    set(gca,'CameraPosition',[-52.1999   30.4749   62.2119])
    set(gca,'CameraUpVector',[1.8643e-001  9.7433e-001 -1.2615e-001])
    set(gca,'CameraViewAngle',6.3060e+000)
    grid on
    set(gcf,'Color','white');
    hl=light('Position',[0,-1,0]);
    set(gca,'Color','white');
    set(gca,'DataAspectRatio',[1,1,1])
    set(gca,'XLim',XLim,'YLim',YLim,'ZLim',ZLim)
    
    HIn=[];
    PhasePlot=linspace(0,20*pi,1500);
    for i=1:length(PhasePlot)
       H=PlotTurbineGeom(T,hf,PhasePlot(i),HIn,Trans,PlotVec,SFVec);
       HIn=H;
       
    end
    
    
    
    
    
    
    
end


