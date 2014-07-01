function H=PlotTurbineGeom(T,HF,Phase,HIn,varargin)

% Plots turbine geometry in structure T to figure HF at rotation phase
% angle Phase (rad). 
% HIn is a cell array containing vectors of handles for any objects to
% update (created on first call and returned in H).
% Optional args:
%   Transparency: Surface transparency value (0 to 1, 1 is opaque, default: .5)
%   PlotVec: 1 to plot element normal and tangent vectors
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
    tx=TR.B(i).tx;
    ty=TR.B(i).ty;
    tz=TR.B(i).tz;
    CtoR=TR.B(i).CtoR;
    
    PEx=TR.B(i).PEx;
    PEy=TR.B(i).PEy;
    PEz=TR.B(i).PEz;
    tEx=TR.B(i).tEx;
    tEy=TR.B(i).tEy;
    tEz=TR.B(i).tEz;
    nEx=TR.B(i).nEx;
    nEy=TR.B(i).nEy;
    nEz=TR.B(i).nEz;

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
            set(HIn{1}{i}(2),'XData',PEx);
            set(HIn{1}{i}(2),'YData',PEy);
            set(HIn{1}{i}(2),'ZData',PEz);
            set(HIn{1}{i}(2),'UData',nEx);
            set(HIn{1}{i}(2),'VData',nEy);
            set(HIn{1}{i}(2),'WData',nEz);
            set(HIn{1}{i}(2),'AutoScaleFactor',SFVec);
            set(HIn{1}{i}(3),'XData',PEx);
            set(HIn{1}{i}(3),'YData',PEy);
            set(HIn{1}{i}(3),'ZData',PEz);
            set(HIn{1}{i}(3),'UData',tEx);
            set(HIn{1}{i}(3),'VData',tEy);
            set(HIn{1}{i}(3),'WData',tEz);
            set(HIn{1}{i}(3),'AutoScaleFactor',SFVec);
        end
        
        drawnow;
        
        H{1}{i}=HIn{1}{i};
    else
        figure(HF)
        hold on
        hs=surf(X,Y,Z);
        set(hs,'FaceColor','g');
        set(hs,'FaceAlpha',SurfTrans);
        set(hs,'LineStyle','-');
        set(hs,'FaceLighting','none')
        set(hs,'EdgeColor',[0 0 0])
        
        % plot normals
        hnv=[];
        htv=[];
        if PlotVec
            hnv=quiver3(PEx,PEy,PEz,nEx,nEy,nEz,SFVec);
            set(hnv,'Color','b')
            htv=quiver3(PEx,PEy,PEz,tEx,tEy,tEz,SFVec);
            set(htv,'Color','g')
        end
        
        % create array of handles for this blade
        H{1}{i}=[hs,hnv,htv];
    end
end

% Surf plot strut geometry
for i=1:TR.NStrut
    MCx=TR.S(i).MCx;
    MCy=TR.S(i).MCy;
    MCz=TR.S(i).MCz;
    CtoR=TR.S(i).CtoR;
    
    % Assume t normal to turbine rotation axis and strut spanwise vector.
    % It's generally assumed that the strut lies in the plane normal to rotation axis,
    % otherwise, should probably be modeled as a blade...
    NElem=TR.S(i).NElem;
    sEx=TR.S(i).sEx;
    sEy=TR.S(i).sEy;
    sEz=TR.S(i).sEz;
    txS=zeros(1,NElem+1);
    tyS=zeros(1,NElem+1);
    tzS=zeros(1,NElem+1);
    for j=1:NElem
        v=[sEx(j),sEy(j),sEz(j)];
        t=cross(TR.RotN,v);
        tMag=sqrt(sum(t.^2));
        txS(j)=t(1)/tMag;
        tyS(j)=t(2)/tMag;
        tzS(j)=t(3)/tMag;
    end
    txS(NElem+1)=txS(NElem);
    tyS(NElem+1)=tyS(NElem);
    tzS(NElem+1)=tzS(NElem);
    
    % Leading and trailing edges
    LEx=MCx-1/2*CtoR.*txS;
    LEy=MCy-1/2*CtoR.*tyS;
    LEz=MCz-1/2*CtoR.*tzS;
    TEx=MCx+1/2*CtoR.*txS;
    TEy=MCy+1/2*CtoR.*tyS;
    TEz=MCz+1/2*CtoR.*tzS;

    % plot strut
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

