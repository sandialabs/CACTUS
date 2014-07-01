function BR=RotateBlade(B,Theta,nR,Origin)

% Rotates blade structure around normal vector nR (size 1 x 3) through
% angle Theta (rad), using specified origin point (size 1 x 3).

% Rotate element end geometry and recalculate element geom for consistency.

% Rotate QC locations
QC=[B.QCx;B.QCy;B.QCz];
QCR=QuatRot(QC',Theta,nR,Origin);
B.QCx=QCR(:,1)';
B.QCy=QCR(:,2)';
B.QCz=QCR(:,3)';

% Rotate t vectors
t=[B.tx;B.ty;B.tz];
tR=QuatRot(t',Theta,nR,[0,0,0]);
B.tx=tR(:,1)';
B.ty=tR(:,2)';
B.tz=tR(:,3)';

% Calc element geometry
B=CalcBEGeom(B);

BR=B;