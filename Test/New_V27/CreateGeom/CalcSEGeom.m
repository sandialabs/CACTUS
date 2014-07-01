function SU=CalcSEGeom(S)
% Calculates the strut element geometry for blade S based on the element end
% geometry in S filled out by the user.
% S: Strut structure to operate on.
%
% Element end geometry fields that must be filled before running this function:
%   NElem
%   MCx
%   MCy
%   MCz
%   CtoR
%
% Element geometry fields calculated by this function:
%   PEx
%   PEy
%   PEz
%   sEx
%   sEy
%   sEz
%   ECtoR
%   EAreaR

for i=1:S.NElem
    PE=[S.MCx(i+1)+S.MCx(i),S.MCy(i+1)+S.MCy(i),S.MCz(i+1)+S.MCz(i)]/2;
    sE=[S.MCx(i+1)-S.MCx(i),S.MCy(i+1)-S.MCy(i),S.MCz(i+1)-S.MCz(i)];
    sEM=sqrt(sum(sE.^2));
    sE=sE./sEM;
    S.PEx(i)=PE(1);
    S.PEy(i)=PE(2);
    S.PEz(i)=PE(3);
    S.sEx(i)=sE(1);
    S.sEy(i)=sE(2);
    S.sEz(i)=sE(3);
    
    % Calc element area and chord
    S.ECtoR(i)=(S.CtoR(i)+S.CtoR(i+1))/2;
    S.EAreaR(i)=sEM*S.ECtoR(i);
end

SU=S;