function T=CreateTurbine(NBlade,NBElem,NStrut,NSElem,RefR,RotN,RotP,RefAR,Type,varargin)

% Creates a CACTUS turbine geometry structure.
%
% NBlade: Number of blades. 
% NBElem: Number of elements per blade.
% NStrut: Number of struts.
% NSElem: Number of elements per strut.
% RefR: Reference radius (ft) defining the scale of the turbine geometry.
%   Will be used to normalize other distances input directly to CACTUS, and
%   also used when outputing dimensional results from CACTUS.
% RotN: Normal vector (size 1 x 3) of turbine rotation axis. Input value
%   used as default when Type is empty, but will be overwritten if a Type 
%   is selected...
% RotP: Origin point (size 1 x 3) on turbine rotation axis. Input value
%   used as default when Type is empty, but will be overwritten if a Type 
%   is selected...
% RefAR: Reference frontal area scaled by RefR^2. Used for
%   force/torque/power coefficient normalization in CACTUS. Input value
%   used as default when Type is empty, but will be overwritten if a Type 
%   is selected...
% Type: Recognized generic turbine type string (see code below). Input empty to 
%   just create an empty turbine structure. If input string is
%   recognized, will fill arrays with actual data for given turbine type
%   using additional arguments defined below.
% varargin: Additional args (comma separated) used for recognized turbine
%   types (see comments in code below for defintion).


T.NBlade=NBlade;
T.NStrut=NStrut;
T.RotN=RotN/sqrt(sum(RotN.^2)); % force normalize
T.RotP=RotP;
T.RefAR=RefAR;
T.RefR=RefR;
T.Type=Type;

% Create blades
for i=1:NBlade
   T.B(i)=CreateBlade(NBElem); 
end

% Create struts
for i=1:NStrut
   T.S(i)=CreateStrut(NSElem); 
end

% Fill geometry if turbine type recognized
if strcmp(Type,'VAWT')==1
    
	% Cross-flow turbine generator for a vertical axis wind turbine (VAWT) with either straight or parabolic blades.
    % For NStrut > 0, struts will be evenly distributed amongst blades (NStrut must be a multiple of NBlade) and 
    % along rotation axis from the center to the tips.
	% Additional arguments:
	% REqR: Equitorial radius to reference radius ratio
	% CR: Blade chord to equitorial radius ratio
	% HR: Turbine height to equitorial radius ratio
	% eta: Blade mount point ratio ((distance behind leading edge of the blade mount point) / (chord))
	% BShape: 0 for straight blades, 1 for parabolic blades
    % CRs: Strut chord to equitorial radius ratio (only used if NStrut > 0)
    % TCs: Strut thickness to chord ratio (only used if NStrut > 0)

    % Get vars
    if length(varargin)<7
        error('Not enough inputs for selected turbine type');
    end 
    REqR=varargin{1};
    CR=varargin{2};
    HR=varargin{3};
    eta=varargin{4};
    BShape=varargin{5};
    CRs=varargin{6};
    TCs=varargin{7};
    
    % Ref to reference radius
    CR=CR*REqR;
    HR=HR*REqR;
    
    % Set rotation axis along y
    T.RotN=[0,1,0];
    T.RotP=[0,0,0];
    
    % Radius ratio function
    yB=linspace(0,HR,NBElem+1);
    if BShape
        % parabolic blades
        rr=REqR*(1-4*(yB/HR-.5).^2);
        % Frontal area normalized by RefR^2
        T.RefAR=2*(REqR*HR-1/3*HR);
        
        % Fill first blade
        deltac=(eta-.25)*CR;
        T.B(1).CtoR=CR*ones(1,NBElem+1);
        T.B(1).tx=ones(1,NBElem+1);
        T.B(1).ty=zeros(1,NBElem+1);
        T.B(1).tz=zeros(1,NBElem+1);
        % normal vector (machine inward)
        drdy=-8/HR*(yB/HR-1/2);
        n=[zeros(1,NBElem+1);drdy;ones(1,NBElem+1)];
        nmag=sqrt(sum(n.^2));
        n=n./nmag(ones(1,3),:);
        Ind=find(n(3,:)<0);
        n(:,Ind)=-n(:,Ind);
        T.B(1).nx=n(1,:);
        T.B(1).ny=n(2,:);
        T.B(1).nz=n(3,:);
        T.B(1).QCx=-deltac*ones(1,NBElem+1);
        T.B(1).QCy=yB;
        T.B(1).QCz=-rr;
    else
        % straight blades
        rr=REqR*ones(size(yB));
        % Frontal area normalized by RefR^2
        T.RefAR=2*REqR*HR;
        
        % Fill first blade
        deltac=(eta-.25)*CR;
        T.B(1).CtoR=CR*ones(1,NBElem+1);
        T.B(1).tx=ones(1,NBElem+1);
        T.B(1).ty=zeros(1,NBElem+1);
        T.B(1).tz=zeros(1,NBElem+1);
        % normal vector (machine inward)
        T.B(1).nx=zeros(1,NBElem+1);
        T.B(1).ny=zeros(1,NBElem+1);
        T.B(1).nz=ones(1,NBElem+1);
        T.B(1).QCx=-deltac*ones(1,NBElem+1);
        T.B(1).QCy=yB;
        T.B(1).QCz=-rr;
    end

    % Calc element area
    T.B(1)=UpdateBElemArea(T.B(1));
    
    % Copy and rotate for other blades
    Phase=linspace(0,2*pi,NBlade+1);
    for i=2:NBlade
        T.B(i)=RotateBlade(T.B(1),Phase(i),T.RotN,T.RotP);
    end
    
    % Fill struts on first blade
    if mod(NStrut,NBlade)~=0
        error('Number of struts must be a multiple of the number of blades for the ''VAWT'' input type.');
    end
    NSpB=round(NStrut/NBlade);
    yS=linspace(0,HR,NSpB+2);
    yS=yS(2:end-1);
    rrS=interp1(yB,rr,yS);
    yC=(yB(2:end)+yB(1:end-1))/2;
    
    for i=1:NSpB
        T.S(i).SEx=zeros(1,NSElem+1);
        T.S(i).SEy=yS(i)*ones(1,NSElem+1);
        T.S(i).SEz=-linspace(0,rrS(i),NSElem+1);
        T.S(i).CtoR=CRs*ones(1,NSElem+1);
        T.S(i).TtoC=TCs;
        T.S(i).BInd=1;
        [m,T.S(i).EInd]=min(abs(yC-yS(i)));
        
        % Calc element area
        T.S(i)=UpdateSElemArea(T.S(i));
    end
    
    % Copy and rotate for other blades
    for i=2:NBlade
        for j=1:NSpB
            SInd=(i-1)*NSpB+j;
            T.S(SInd)=RotateStrut(T.S(j),Phase(i),T.RotN,T.RotP);
            T.S(i).BInd=i;
        end
    end

    
elseif strcmp(Type,'HAWT')==1
    
    % Axial-flow turbine generator for a horizontal axis wind turbine (HAWT).
	% Additional arguments:
    % RMaxR: Turbine radius to reference radius ratio
    % HubRR: Hub radius to turbine radius ratio
    % CR: Blade chord to turbine radius ratio (NBElem+1 elements ordered root to tip)
    % bTwist: Blade planform twist at each element end (deg, w.r.t. rotor disk plane, positive LE into the wind (-x), NBElem+1 elements ordered root to tip)
    % bi: Blade planform incidence (deg, w.r.t. rotor disk plane, positive LE into the wind (-x))
    % eta: Blade mount point ratio ((distance behind leading edge of the blade mount point) / (chord))
    % bCone: Blade coning angle (deg, positive tip into the wind (-x))
    % Tilt: Rotor tilt angle (deg, positive windward axis tilted up)

    % Get vars
    if length(varargin)<8
        error('Not enough inputs for selected turbine type');
    end 
    RMaxR=varargin{1};
    HubRR=varargin{2};
    CR=varargin{3};
    bTwist=varargin{4};
    bi=varargin{5};
    eta=varargin{6};
    bCone=varargin{7};
    Tilt=varargin{8};
    
    % Ref to reference radius
    CR=CR*RMaxR;
    HubRR=HubRR*RMaxR;
    
    % Set rotation axis along x
    T.RotN=[1,0,0];
    T.RotP=[0,0,0];
    
    % Radius ratio function
    rB=linspace(HubRR,RMaxR,NBElem+1);
    % Frontal area normalized by RefR^2
    T.RefAR=pi*RMaxR^2;

    % Fill first blade
    deltac=(eta-.25)*CR(1);
    T.B(1).QCx=zeros(1,NBElem+1);
    T.B(1).QCy=rB;
    T.B(1).QCz=deltac*ones(1,NBElem+1);
    T.B(1).CtoR=CR;
    sTwist=sin(bTwist/180*pi);
    cTwist=cos(bTwist/180*pi);
    T.B(1).tx=sTwist;
    T.B(1).ty=zeros(1,NBElem+1);
    T.B(1).tz=-cTwist;
    % normal vector (machine rearward (x))
    T.B(1).nx=cTwist;
    T.B(1).ny=zeros(1,NBElem+1);
    T.B(1).nz=sTwist;

    % Calc element area
    T.B(1)=UpdateBElemArea(T.B(1));
    
    % Rotate through incidence and coning angle
    T.B(1)=RotateBlade(T.B(1),bi/180*pi,[0,-1,0],[0,0,0]);
    T.B(1)=RotateBlade(T.B(1),bCone/180*pi,[0,0,1],[0,0,0]);
    
    % Copy and rotate for other blades
    Phase=linspace(0,2*pi,NBlade+1);
    for i=2:NBlade
        T.B(i)=RotateBlade(T.B(1),Phase(i),T.RotN,T.RotP);
    end
    
    % Rotate turbine through tilt angle
    T=RotateTurbine(T,Tilt/180*pi,[0,0,-1],[0,0,0]);
    
end
