function BU=CalcBEGeom(B)
% Calculates the blade element geometry for blade B based on the element end
% geometry in B filled out by the user.
% B: Blade structure to operate on.
%
% Element end geometry fields that must be filled before running this function:
%   NElem
%   FlipN
%   QCx
%   QCy
%   QCz
%   tx
%   ty
%   tz
%   CtoR
%
% Element geometry fields calculated by this function:
%   PEx
%   PEy
%   PEz
%   tEx
%   tEy
%   tEz
%   nEx
%   nEy
%   nEz
%   sEx
%   sEy
%   sEz
%   ECtoR
%   EAreaR

NElem=B.NElem;
FlipN=B.FlipN;

for i=1:NElem
    PE=[B.QCx(i+1)+B.QCx(i),B.QCy(i+1)+B.QCy(i),B.QCz(i+1)+B.QCz(i)]/2;
    sE=-[B.QCx(i+1)-B.QCx(i),B.QCy(i+1)-B.QCy(i),B.QCz(i+1)-B.QCz(i)]; % nominal element spanwise direction set opposite to QC line in CACTUS
    sEM=sqrt(sum(sE.^2));
    sE=sE./sEM;
    tE=[B.tx(i+1)+B.tx(i),B.ty(i+1)+B.ty(i),B.tz(i+1)+B.tz(i)]/2;
    % Force tE normal to sE
    tE=tE-(tE*sE')*sE;
    tEM=sqrt(sum(tE.^2));
    if tEM<1e-10
        error('Error: Element t vector must not be parallel to quarter chord line.')
    end
    tE=tE./tEM;
    % Calc normal vector
    nE=cross(sE,tE);
    nE=nE./sqrt(sum(nE.^2));
    
    % Flip normal direction if requested
    if FlipN==1
        nE=-nE;
        sE=-sE;
    end
    
    B.PEx(i)=PE(1);
    B.PEy(i)=PE(2);
    B.PEz(i)=PE(3);
    B.tEx(i)=tE(1);
    B.tEy(i)=tE(2);
    B.tEz(i)=tE(3);
    B.nEx(i)=nE(1);
    B.nEy(i)=nE(2);
    B.nEz(i)=nE(3);
    B.sEx(i)=sE(1);
    B.sEy(i)=sE(2);
    B.sEz(i)=sE(3);
    
    % Calc element area and chord
    S=[-1/4,3/4];
    SR=[3/4,-1/4];
    % Calc quad area from two triangular facets
    Px=[B.QCx(i)+S*B.CtoR(i)*B.tx(i),B.QCx(i+1)+SR*B.CtoR(i+1)*B.tx(i+1),B.QCx(i)-1/4*B.CtoR(i)*B.tx(i)];
    Py=[B.QCy(i)+S*B.CtoR(i)*B.ty(i),B.QCy(i+1)+SR*B.CtoR(i+1)*B.ty(i+1),B.QCy(i)-1/4*B.CtoR(i)*B.ty(i)];
    Pz=[B.QCz(i)+S*B.CtoR(i)*B.tz(i),B.QCz(i+1)+SR*B.CtoR(i+1)*B.tz(i+1),B.QCz(i)-1/4*B.CtoR(i)*B.tz(i)];
    V=[diff(Px);diff(Py);diff(Pz)];
    A1=cross(V(:,1),V(:,2))/2;
    A2=cross(V(:,3),V(:,4))/2;
    B.EAreaR(i)=sqrt(sum(A1.^2))+sqrt(sum(A2.^2));
    % Calc average element chord from area and span
    B.ECtoR(i)=B.EAreaR(i)/sEM;
end

BU=B;