function B=CreateBlade(NElem)

% Create empty blade with specified number of elements, NElem. Quarter
% chord line and section normal and tangential vectors defined at zero
% turbine rotation phase.

% QC: Quarter chord location (over ref. radius) for each element end (size NElem + 1).
% n: Blade section normal vector for each element end (size NElem + 1). Defines
%   positive angle of attack (AOA positive when relative velocity component
%   is positive in the normal direction.
% t: Blade section tangent vector (must use rearward chord line) for each element end
%   (size NElem + 1).
% CtoR: Chord to ref. radius for each element end (size NElem + 1).
% AreaR: Element area over ref. radius squared for each element (size NElem).
% iSect: Airfoil section index for each element (size NElem). Used in
%   CACTUS to identify the airfoil data tables (defined in the CACTUS input
%   file) to use with that element.

B.NElem=NElem;
B.QCx=zeros(1,NElem+1);
B.QCy=zeros(1,NElem+1);
B.QCz=zeros(1,NElem+1);
B.nx=zeros(1,NElem+1);
B.ny=zeros(1,NElem+1);
B.nz=zeros(1,NElem+1);
B.tx=zeros(1,NElem+1);
B.ty=zeros(1,NElem+1);
B.tz=zeros(1,NElem+1);
B.CtoR=zeros(1,NElem+1);
B.AreaR=zeros(1,NElem);
B.iSect=ones(1,NElem);
