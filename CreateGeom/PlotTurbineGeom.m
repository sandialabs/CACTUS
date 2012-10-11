function H=PlotTurbineGeom(T,HF,Phase,HIn,varargin)

% Plots turbine geometry in structure T to figure HF at rotation phase
% angle Phase (rad). 
% HIn is a cell array containing vectors of handles for any objects to
% update (created on first call and returned in H).
% Optional args:
%   PlotVec: 1 to plot element normal vectors
%   SFVec: scale factor for normal vector plotting (default: 1)

PlotVec=0;
SFVec=1;
if nargin==5
    PlotVec=varargin{1};
elseif nargin==6
    PlotVec=varargin{1};
    SFVec=varargin{2};
end

% Rotate turbine to phase
nR=T.RotN;
oR=T.RotP;
TR=RotateTurbine(T,Phase,nR,oR);

% Surf plot blade geometry
for i=1:length(TR.B)
    QCx=TR.B(i).QCx;
    QCy=TR.B(i).QCy;
    QCz=TR.B(i).QCz;
    nx=TR.B(i).nx;
    ny=TR.B(i).ny;
    nz=TR.B(i).nz;
    tx=TR.B(i).tx;
    ty=TR.B(i).ty;
    tz=TR.B(i).tz;
    CtoR=TR.B(i).CtoR;
    
    % Leading and trailing edges
    LEx=QCx-1/4*CtoR.*tx;
    LEy=QCy-1/4*CtoR.*ty;
    LEz=QCz-1/4*CtoR.*tz;
    TEx=QCx+3/4*CtoR.*tx;
    TEy=QCy+3/4*CtoR.*ty;
    TEz=QCz+3/4*CtoR.*tz;

    % plot blade
    X=[LEx;TEx];
    Y=[LEy;TEy];
    Z=[LEz;TEz];
    
    % If plot objects exist, reset data and redraw for smoother
    % animation. Otherwise plot new...
    if ~isempty(HIn)
        set(HIn{i}(1),'XData',X);
        set(HIn{i}(1),'YData',Y);
        set(HIn{i}(1),'ZData',Z);
        
        if PlotVec
            set(HIn{i}(2),'XData',QCx);
            set(HIn{i}(2),'YData',QCy);
            set(HIn{i}(2),'ZData',QCz);
            set(HIn{i}(2),'UData',nx);
            set(HIn{i}(2),'VData',ny);
            set(HIn{i}(2),'WData',nz);
            set(HIn{i}(2),'AutoScaleFactor',SFVec);
        end
        
        drawnow;
        
        H{i}=HIn{i};
    else
        figure(HF)
        hold on
        hs=surf(X,Y,Z);
        set(hs,'FaceColor','b');
        set(hs,'FaceAlpha',.5);
        set(hs,'LineStyle','none');
        set(hs,'FaceLighting','none')
        
        % plot normals
        hv=[];
        if PlotVec
            hv=quiver3(QCx,QCy,QCz,nx,ny,nz,SFVec);
            set(hv,'Color','b')
        end
        
        % create array of handles for this blade
        H{i}=[hs,hv];
    end
end

