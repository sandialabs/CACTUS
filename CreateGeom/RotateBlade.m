function BR=RotateBlade(B,Theta,nR,Origin)

% Rotates blade structure around normal vector nR (size 1 x 3) through
% angle Theta (rad), using specified origin point (size 1 x 3).

% Rotate QC locations
QC=[B.QCx;B.QCy;B.QCz];
QCR=QuatRot(QC',Theta,nR,Origin);
B.QCx=QCR(:,1)';
B.QCy=QCR(:,2)';
B.QCz=QCR(:,3)';

% Rotate n and t vectors
t=[B.tx;B.ty;B.tz];
tR=QuatRot(t',Theta,nR,[0,0,0]);
B.tx=tR(:,1)';
B.ty=tR(:,2)';
B.tz=tR(:,3)';

n=[B.nx;B.ny;B.nz];
nR=QuatRot(n',Theta,nR,[0,0,0]);
B.nx=nR(:,1)';
B.ny=nR(:,2)';
B.nz=nR(:,3)';

BR=B;