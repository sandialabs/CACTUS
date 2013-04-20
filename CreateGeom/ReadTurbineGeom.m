function T=ReadTurbineGeom(FN)

% Read a CACTUS geometry input file.

fid=fopen(FN,'r');

S=fgetl(fid);
Ind=strfind(S,':');
A=sscanf(S(Ind+1:end),'%d');
T.NBlade=A;

S=fgetl(fid);
Ind=strfind(S,':');
A=sscanf(S(Ind+1:end),'%d');
T.NStrut=A;

S=fgetl(fid);
Ind=strfind(S,':');
A=sscanf(S(Ind+1:end),'%e');
T.RotN=A';

S=fgetl(fid);
Ind=strfind(S,':');
A=sscanf(S(Ind+1:end),'%e');
T.RotP=A';

S=fgetl(fid);
Ind=strfind(S,':');
A=sscanf(S(Ind+1:end),'%e');
T.RefAR=A;

S=fgetl(fid);
Ind=strfind(S,':');
A=sscanf(S(Ind+1:end),'%e');
T.RefR=A;

S=fgetl(fid);
Ind=strfind(S,':');
T.Type=strtrim(S(Ind+1:end));

for i=1:T.NBlade
    S=fgetl(fid);
    
    S=fgetl(fid);
    Ind=strfind(S,':');
    A=sscanf(S(Ind+1:end),'%d');
    T.B(i).NElem=A;
    
    S=fgetl(fid);
    Ind=strfind(S,':');
    A=sscanf(S(Ind+1:end),'%e');
    T.B(i).QCx=A';

    S=fgetl(fid);
    Ind=strfind(S,':');
    A=sscanf(S(Ind+1:end),'%e');
    T.B(i).QCy=A';
    
    S=fgetl(fid);
    Ind=strfind(S,':');
    A=sscanf(S(Ind+1:end),'%e');
    T.B(i).QCz=A';
    
    S=fgetl(fid);
    Ind=strfind(S,':');
    A=sscanf(S(Ind+1:end),'%e');
    T.B(i).nx=A';

    S=fgetl(fid);
    Ind=strfind(S,':');
    A=sscanf(S(Ind+1:end),'%e');
    T.B(i).ny=A';    

    S=fgetl(fid);
    Ind=strfind(S,':');
    A=sscanf(S(Ind+1:end),'%e');
    T.B(i).nz=A';

    S=fgetl(fid);
    Ind=strfind(S,':');
    A=sscanf(S(Ind+1:end),'%e');
    T.B(i).tx=A';

    S=fgetl(fid);
    Ind=strfind(S,':');
    A=sscanf(S(Ind+1:end),'%e');
    T.B(i).ty=A';    

    S=fgetl(fid);
    Ind=strfind(S,':');
    A=sscanf(S(Ind+1:end),'%e');
    T.B(i).tz=A';
    
    S=fgetl(fid);
    Ind=strfind(S,':');
    A=sscanf(S(Ind+1:end),'%e');
    T.B(i).CtoR=A';

    S=fgetl(fid);
    Ind=strfind(S,':');
    A=sscanf(S(Ind+1:end),'%e');
    T.B(i).AreaR=A';
    
    S=fgetl(fid);
    Ind=strfind(S,':');
    A=sscanf(S(Ind+1:end),'%d');
    T.B(i).iSect=A';
end

for i=1:T.NStrut
    S=fgetl(fid);
    
    S=fgetl(fid);
    Ind=strfind(S,':');
    A=sscanf(S(Ind+1:end),'%d');
    T.S(i).NElem=A;
    
    S=fgetl(fid);
    Ind=strfind(S,':');
    A=sscanf(S(Ind+1:end),'%e');
    T.S(i).SEx=A';

    S=fgetl(fid);
    Ind=strfind(S,':');
    A=sscanf(S(Ind+1:end),'%e');
    T.S(i).SEy=A';
    
    S=fgetl(fid);
    Ind=strfind(S,':');
    A=sscanf(S(Ind+1:end),'%e');
    T.S(i).SEz=A';
    
    S=fgetl(fid);
    Ind=strfind(S,':');
    A=sscanf(S(Ind+1:end),'%e');
    T.S(i).CtoR=A';

    S=fgetl(fid);
    Ind=strfind(S,':');
    A=sscanf(S(Ind+1:end),'%e');
    T.S(i).AreaR=A';
    
    S=fgetl(fid);
    Ind=strfind(S,':');
    A=sscanf(S(Ind+1:end),'%e');
    T.S(i).TtoC=A;
    
    S=fgetl(fid);
    Ind=strfind(S,':');
    A=sscanf(S(Ind+1:end),'%d');
    T.S(i).BIndS=A;
    
    S=fgetl(fid);
    Ind=strfind(S,':');
    A=sscanf(S(Ind+1:end),'%d');
    T.S(i).EIndS=A;
    
    S=fgetl(fid);
    Ind=strfind(S,':');
    A=sscanf(S(Ind+1:end),'%d');
    T.S(i).BIndE=A;
    
    S=fgetl(fid);
    Ind=strfind(S,':');
    A=sscanf(S(Ind+1:end),'%d');
    T.S(i).EIndE=A;    
end

fclose(fid);

% Force normalize unit vectors
T.RotN=T.RotN/sqrt(sum(T.RotN.^2));
for i=1:T.NBlade
    % normal
    nv=[T.B(i).nx;T.B(i).ny;T.B(i).nz];
    nvMag=sqrt(sum(nv.^2));
    nv=nv./nvMag(ones(3,1),:);
    T.B(i).nx=nv(1,:);
    T.B(i).ny=nv(2,:);
    T.B(i).nz=nv(3,:);
    
    % tangential
    tv=[T.B(i).tx;T.B(i).ty;T.B(i).tz];
    tvMag=sqrt(sum(tv.^2));
    tv=tv./tvMag(ones(3,1),:);
    T.B(i).tx=tv(1,:);
    T.B(i).ty=tv(2,:);
    T.B(i).tz=tv(3,:);
end
