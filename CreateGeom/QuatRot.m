function vR=QuatRot(v,Theta,nR,Origin)

% Perform rotation of vector(s) v around normal vector nR using the
% quaternion machinery.
% v: input vector(s) (can be m x 3 for m vectors to rotate)
% Theta: rotation angle (rad)
% nR: normal vector around which to rotate
% Origin: origin point of rotation
%
% vR: Rotated vector(s) (size m x 3 for m input vectors)

% Force normalize nR
nR=nR/sqrt(sum(nR.^2));

% Quaternion form of v
O=Origin(ones(size(v,1),1),:);
vO=v-O;
p=[zeros(size(v,1),1),vO];

% Rotation quaternion and conjugate
q=[cos(Theta/2),nR*sin(Theta/2)];
qbar=[q(1),-q(2:4)];

QL=[q(1) -q(2) -q(3) -q(4)
    q(2)  q(1) -q(4)  q(3)
    q(3)  q(4)  q(1) -q(2)
    q(4) -q(3)  q(2)  q(1)];
   
QbarR=[qbar(1) -qbar(2) -qbar(3) -qbar(4)
       qbar(2)  qbar(1)  qbar(4) -qbar(3)
       qbar(3) -qbar(4)  qbar(1)  qbar(2)
       qbar(4)  qbar(3) -qbar(2)  qbar(1)];

% Rotate p
pR=p*(QbarR*QL)';
vR=pR(:,2:4)+O;



