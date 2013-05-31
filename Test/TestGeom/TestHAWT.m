clear
close all

% Creates test HAWT geometry file

% Add geom creation scripts to path
path(path,'/ascldap/users/jmurray/Project/CACTUS/stable/trunk/CreateGeom');

% Params
R=18.15; % radius (ft)
HubRR=0.0334; % hub radius ratio
Tilt=0; % Rotor tilt angle (deg, positive windward axis tilted up)
eta=.3; % Blade mount point ratio ((distance behind leading edge of the blade mount point) / (root chord))
CRr=[0.0394,0.1332,0.1133,0.0926,0.0759,0.0551]; % chord to radius
bCone=0; % Blade coning angle (deg, positive tip into the wind (-x))
bi=8.5; % Blade planform incidence (deg, w.r.t. rotor disk plane, positive LE into the wind (-x))
bTwist=[0,20,4.715,0.494,-0.92,-2.5]; % Blade planform twist at each element end (same sign as bi, elements ordered root to tip)	
NBlade=2;
NElem=5;

% Output filename
FN='TestHAWT.geom';

% Plot data?
PlotTurbine=1;

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
    XLim=[-2,2];
    YLim=[-2,2];
    ZLim=[-2,2];
    
    % Plot element normals
    PlotVec=1;
    SFVec=.5;
    
    hf=figure(1);
    set(hf,'Position',[303   124   956   610]) 
    set(gca,'Position',[5.2743e-002  5.1245e-002  8.9979e-001  8.8141e-001])
    set(gca,'CameraPosition',[-52.1999   30.4749   62.2119])
    set(gca,'CameraUpVector',[1.8643e-001  9.7433e-001 -1.2615e-001])
    set(gca,'CameraViewAngle',6.3060e+000)
    grid on
    set(gcf,'Color','white');
    hl=light('Position',[-1,0,0]);
    set(gca,'Color','white');
    set(gca,'DataAspectRatio',[1,1,1])
    set(gca,'XLim',XLim,'YLim',YLim,'ZLim',ZLim)
    
    HIn=[];
    PhasePlot=linspace(0,2*pi,150);
    for i=1:length(PhasePlot)
       H=PlotTurbineGeom(T,hf,PhasePlot(i),HIn,PlotVec,SFVec);
       HIn=H;
       pause(.01);
    end
    
end


