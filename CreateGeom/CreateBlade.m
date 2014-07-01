function B=CreateBlade(NElem)
% Create empty blade with specified number of elements, NElem. Quarter
% chord line and section normal and tangential vectors defined at zero
% turbine rotation phase.
%
% NElem: Number of elements
% FlipN: Flip element normal vector direction flag. A value of 1 indicates flipped
%   direction, 0 indicates nominal direction. Nominal element normal vector direction
%   defined by t vector and the quarter chord line vector (positive in the
%   direction of increasing element number) as n=cross(t,QCv). Note that AOA
%   is defined positive on an element when relative velocity is positive in
%   the normal direction.
% QC: Quarter chord location (over ref. radius) for each element end (size NElem + 1).
% t: Blade section tangent vector (must use rearward chord line) for each element end
%   (size NElem + 1).
% CtoR: Chord to ref. radius for each element end (size NElem + 1).
% PE: Location (over ref. radius) for each element center (size NElem).
% tE: Element tangent vector (size NElem).
% nE: Elememt normal vector (size NElem).
% sE: Element spanwise vector (size NElem).
% ECtoR: Chord to ref. radius for each element (size NElem).
% EAreaR: Element area over ref. radius squared for each element (size NElem).
% iSect: Airfoil section index for each element (size NElem). Used in
%   CACTUS to identify the airfoil data tables (defined in the CACTUS input
%   file) to use with that element.

B.NElem=NElem;
B.FlipN=0;
% Element end geometry
B.QCx=zeros(1,NElem+1);
B.QCy=zeros(1,NElem+1);
B.QCz=zeros(1,NElem+1);
B.tx=zeros(1,NElem+1);
B.ty=zeros(1,NElem+1);
B.tz=zeros(1,NElem+1);
B.CtoR=zeros(1,NElem+1);
% Element geometry
B.PEx=zeros(1,NElem);
B.PEy=zeros(1,NElem);
B.PEz=zeros(1,NElem);
B.tEx=zeros(1,NElem);
B.tEy=zeros(1,NElem);
B.tEz=zeros(1,NElem);
B.nEx=zeros(1,NElem);
B.nEy=zeros(1,NElem);
B.nEz=zeros(1,NElem);
B.sEx=zeros(1,NElem);
B.sEy=zeros(1,NElem);
B.sEz=zeros(1,NElem);
B.ECtoR=zeros(1,NElem);
B.EAreaR=zeros(1,NElem);
B.iSect=ones(1,NElem);
