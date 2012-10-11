function T=CreateTurbine(NBlade,NElem,RefR,RotN,RotP,RefAR,Type,varargin)

% Creates a CACTUS turbine geometry structure.
%
% NBlade: Number of blades. 
% NElem: Number of elements per blade. 
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
% Type: Recognized generic turbine type string (see below). Input empty to 
%   just create an empty turbine structure (default). If input string is
%   recognized, will fill blades with actual data for given turbine type
%   using additional arguments defined below.
% Additional args used for recognized turbine types:
%   Type: 'VAWT'
%       REqR: Equitorial radius to reference radius ratio
%       CR: Blade chord to equitorial radius ratio
%       HR: Turbine height to equitorial radius ratio
%       eta: Blade mount point ratio ((distance behind leading edge of the blade mount point) / (chord))
%       BShape: 0 for straight blades, 1 for parabolic blades
%   Type: 'HAWT'
%       RMaxR: Turbine radius to reference radius ratio
%       HubRR: Hub radius to turbine radius ratio
%       CR: Blade chord to turbine radius ratio (NElem+1 elements ordered root to tip)
%       bTwist: Blade planform twist at each element end (same sign as bi, NElem+1 elements ordered root to tip)
%       eta: Blade mount point ratio ((distance behind leading edge of the blade mount point) / (chord))
%       bCone: Blade coning angle (deg, positive tip into the wind (-x))
%       Tilt: Rotor tilt angle (deg, positive windward axis tilted up)


T.NBlade=NBlade;
T.RotN=RotN/sqrt(sum(RotN.^2)); % force normalize
T.RotP=RotP;
T.RefAR=RefAR;
T.RefR=RefR;
T.Type=Type;

% Create blades
for i=1:NBlade
   T.B(i)=CreateBlade(NElem); 
end

% Fill blades if turbine type recognized
if strcmp(Type,'VAWT')==1
    
	% Cross-flow turbine generator for a vertical axis wind turbine (VAWT) with either straight or parabolic blades.
	% Additional arguments:
	% REqR: Equitorial radius to reference radius ratio
	% CR: Blade chord to equitorial radius ratio
	% HR: Turbine height to equitorial radius ratio
	% eta: Blade mount point ratio ((distance behind leading edge of the blade mount point) / (chord))
	% BShape: 0 for straight blades, 1 for parabolic blades

    % Get vars
    if length(varargin)<5
        error('Not enough inputs for selected turbine type');
    end 
    REqR=varargin{1};
    CR=varargin{2};
    HR=varargin{3};
    eta=varargin{4};
    BShape=varargin{5};
    
    % Ref to reference radius
    CR=CR*REqR;
    HR=HR*REqR;
    
    % Set rotation axis along y
    T.RotN=[0,1,0];
    T.RotP=[0,0,0];
    
    % Radius ratio function
    yB=linspace(0,HR,NElem+1);
    if BShape
        % parabolic blades
        rr=REqR*(1-4*(yB/HR-.5).^2);
        % Frontal area normalized by RefR^2
        T.RefAR=2*(REqR*HR-1/3*HR);
        
        % Fill first blade
        deltac=(eta-.25)*CR;
        T.B(1).CtoR=CR*ones(1,NElem+1);
        T.B(1).tx=ones(1,NElem+1);
        T.B(1).ty=zeros(1,NElem+1);
        T.B(1).tz=zeros(1,NElem+1);
        % normal vector (machine inward)
        drdy=-8/HR*(yB/HR-1/2);
        n=[zeros(1,NElem+1);drdy;ones(1,NElem+1)];
        nmag=sqrt(sum(n.^2));
        n=n./nmag(ones(1,3),:);
        Ind=find(n(3,:)<0);
        n(:,Ind)=-n(:,Ind);
        T.B(1).nx=n(1,:);
        T.B(1).ny=n(2,:);
        T.B(1).nz=n(3,:);
        T.B(1).QCx=-deltac*ones(1,NElem+1);
        T.B(1).QCy=yB;
        T.B(1).QCz=-rr;
    else
        % straight blades
        rr=REqR*ones(size(yB));
        % Frontal area normalized by RefR^2
        T.RefAR=2*REqR*HR;
        
        % Fill first blade
        deltac=(eta-.25)*CR;
        T.B(1).CtoR=CR*ones(1,NElem+1);
        T.B(1).tx=ones(1,NElem+1);
        T.B(1).ty=zeros(1,NElem+1);
        T.B(1).tz=zeros(1,NElem+1);
        % normal vector (machine inward)
        T.B(1).nx=zeros(1,NElem+1);
        T.B(1).ny=zeros(1,NElem+1);
        T.B(1).nz=ones(1,NElem+1);
        T.B(1).QCx=-deltac*ones(1,NElem+1);
        T.B(1).QCy=yB;
        T.B(1).QCz=-rr;
    end

    % Calc element area
    T.B(1)=UpdateElemArea(T.B(1));
    
    % Copy and rotate for other blades
    Phase=linspace(0,2*pi,NBlade+1);
    for i=2:NBlade
        T.B(i)=RotateBlade(T.B(1),Phase(i),T.RotN,T.RotP);
    end

elseif strcmp(Type,'HAWT')==1
    
    % Axial-flow turbine generator for a horizontal axis wind turbine (HAWT).
	% Additional arguments:
    % RMaxR: Turbine radius to reference radius ratio
    % HubRR: Hub radius to turbine radius ratio
    % CR: Blade chord to turbine radius ratio (NElem+1 elements ordered root to tip)
    % bTwist: Blade planform twist at each element end (deg, w.r.t. rotor disk plane, positive LE into the wind (-x), NElem+1 elements ordered root to tip)
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
    rB=linspace(HubRR,RMaxR,NElem+1);
    % Frontal area normalized by RefR^2
    T.RefAR=pi*RMaxR^2;

    % Fill first blade
    deltac=(eta-.25)*CR(1);
    T.B(1).QCx=zeros(1,NElem+1);
    T.B(1).QCy=rB;
    T.B(1).QCz=deltac*ones(1,NElem+1);
    T.B(1).CtoR=CR;
    sTwist=sin(bTwist/180*pi);
    cTwist=cos(bTwist/180*pi);
    T.B(1).tx=sTwist;
    T.B(1).ty=zeros(1,NElem+1);
    T.B(1).tz=-cTwist;
    % normal vector (machine rearward (x))
    T.B(1).nx=cTwist;
    T.B(1).ny=zeros(1,NElem+1);
    T.B(1).nz=sTwist;

    % Calc element area
    T.B(1)=UpdateElemArea(T.B(1));
    
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
