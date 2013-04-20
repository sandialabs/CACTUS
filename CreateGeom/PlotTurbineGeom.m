function H=PlotTurbineGeom(T,HF,Phase,HIn,varargin)

% Plots turbine geometry in structure T to figure HF at rotation phase
% angle Phase (rad). 
% HIn is a cell array containing vectors of handles for any objects to
% update (created on first call and returned in H).
% Optional args:
%   Transparency: Surface transparency value (0 to 1, 1 is opaque, default: .5)
%   PlotVec: 1 to plot element normal vectors
%   SFVec: scale factor for normal vector plotting (default: 1)

SurfTrans=.5;
PlotVec=0;
SFVec=1;
if nargin==5
    SurfTrans=varargin{1};
elseif nargin==6
    SurfTrans=varargin{1};
    PlotVec=varargin{2};
elseif nargin == 7 
    SurfTrans=varargin{1};
    PlotVec=varargin{2};
    SFVec=varargin{3};
end

% Rotate turbine to phase
nR=T.RotN;
oR=T.RotP;
TR=RotateTurbine(T,Phase,nR,oR);

% Surf plot blade geometry
for i=1:TR.NBlade
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
        set(HIn{1}{i}(1),'XData',X);
        set(HIn{1}{i}(1),'YData',Y);
        set(HIn{1}{i}(1),'ZData',Z);
        
        if PlotVec
            set(HIn{1}{i}(2),'XData',QCx);
            set(HIn{1}{i}(2),'YData',QCy);
            set(HIn{1}{i}(2),'ZData',QCz);
            set(HIn{1}{i}(2),'UData',nx);
            set(HIn{1}{i}(2),'VData',ny);
            set(HIn{1}{i}(2),'WData',nz);
            set(HIn{1}{i}(2),'AutoScaleFactor',SFVec);
        end
        
        drawnow;
        
        H{1}{i}=HIn{1}{i};
    else
        figure(HF)
        hold on
        hs=surf(X,Y,Z);
        set(hs,'FaceColor','b');
        set(hs,'FaceAlpha',SurfTrans);
        set(hs,'LineStyle','none');
        set(hs,'FaceLighting','none')
        
        % plot normals
        hv=[];
        if PlotVec
            hv=quiver3(QCx,QCy,QCz,nx,ny,nz,SFVec);
            set(hv,'Color','b')
        end
        
        % create array of handles for this blade
        H{1}{i}=[hs,hv];
    end
end

% Surf plot strut geometry
for i=1:TR.NStrut
    SEx=TR.S(i).SEx;
    SEy=TR.S(i).SEy;
    SEz=TR.S(i).SEz;
    CtoR=TR.S(i).CtoR;
    
    % Assume t normal to turbine rotation axis and strut spanwise vector.
    % It's generally assumed that the strut lies in the plane normal to rotation axis,
    % otherwise, should probably be modeled as a blade...
    txS=zeros(1,TR.S(i).NElem+1);
    tyS=zeros(1,TR.S(i).NElem+1);
    tzS=zeros(1,TR.S(i).NElem+1);
    for j=1:TR.S(i).NElem
        v=[SEx(j+1)-SEx(j),SEy(j+1)-SEy(j),SEz(j+1)-SEz(j)];
        t=cross(TR.RotN,v);
        tMag=sqrt(sum(t.^2));
        txS(j)=t(1)/tMag;
        tyS(j)=t(2)/tMag;
        tzS(j)=t(3)/tMag;
    end
    txS(TR.S(i).NElem+1)=txS(TR.S(i).NElem);
    tyS(TR.S(i).NElem+1)=tyS(TR.S(i).NElem);
    tzS(TR.S(i).NElem+1)=tzS(TR.S(i).NElem);
    
    % Leading and trailing edges
    LEx=SEx-1/2*CtoR.*txS;
    LEy=SEy-1/2*CtoR.*tyS;
    LEz=SEz-1/2*CtoR.*tzS;
    TEx=SEx+1/2*CtoR.*txS;
    TEy=SEy+1/2*CtoR.*tyS;
    TEz=SEz+1/2*CtoR.*tzS;

    % plot blade
    X=[LEx;TEx];
    Y=[LEy;TEy];
    Z=[LEz;TEz];
    
    % If plot objects exist, reset data and redraw for smoother
    % animation. Otherwise plot new...
    if ~isempty(HIn)
        set(HIn{2}{i},'XData',X);
        set(HIn{2}{i},'YData',Y);
        set(HIn{2}{i},'ZData',Z);
        
        drawnow;
        
        H{2}{i}=HIn{2}{i};
    else
        figure(HF)
        hold on
        hs=surf(X,Y,Z);
        set(hs,'FaceColor','k');
        set(hs,'FaceAlpha',SurfTrans);
        set(hs,'LineStyle','none');
        set(hs,'FaceLighting','none')
        
        % create array of handles for this blade
        H{2}{i}=hs;
    end
end

