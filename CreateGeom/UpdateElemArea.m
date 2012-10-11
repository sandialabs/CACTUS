function BU=UpdateElemArea(B)

% Updates AreaR component of given blade structure (element areas). Must
% set quarter chord location, tangential vector, and chord to radius
% ratio components in B before running this function.

S=[-1/4,3/4];
SR=[3/4,-1/4];
for i=1:length(B.AreaR)
    % Calc quad area from two triangular facets
    Px=[B.QCx(i)+S*B.CtoR(i)*B.tx(i),B.QCx(i+1)+SR*B.CtoR(i+1)*B.tx(i+1),B.QCx(i)-1/4*B.CtoR(i)*B.tx(i)];
    Py=[B.QCy(i)+S*B.CtoR(i)*B.ty(i),B.QCy(i+1)+SR*B.CtoR(i+1)*B.ty(i+1),B.QCy(i)-1/4*B.CtoR(i)*B.ty(i)];
    Pz=[B.QCz(i)+S*B.CtoR(i)*B.tz(i),B.QCz(i+1)+SR*B.CtoR(i+1)*B.tz(i+1),B.QCz(i)-1/4*B.CtoR(i)*B.tz(i)];
    V=[diff(Px);diff(Py);diff(Pz)];
    A1=cross(V(:,1),V(:,2))/2;
    A2=cross(V(:,3),V(:,4))/2;
    B.AreaR(i)=sqrt(sum(A1.^2))+sqrt(sum(A2.^2));
end

BU=B;