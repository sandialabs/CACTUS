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
    A=sscanf(S(Ind+1:end),'%d');
    T.B(i).FlipN=A;
    
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
    T.B(i).PEx=A';

    S=fgetl(fid);
    Ind=strfind(S,':');
    A=sscanf(S(Ind+1:end),'%e');
    T.B(i).PEy=A';
    
    S=fgetl(fid);
    Ind=strfind(S,':');
    A=sscanf(S(Ind+1:end),'%e');
    T.B(i).PEz=A';

    S=fgetl(fid);
    Ind=strfind(S,':');
    A=sscanf(S(Ind+1:end),'%e');
    T.B(i).tEx=A';

    S=fgetl(fid);
    Ind=strfind(S,':');
    A=sscanf(S(Ind+1:end),'%e');
    T.B(i).tEy=A';    

    S=fgetl(fid);
    Ind=strfind(S,':');
    A=sscanf(S(Ind+1:end),'%e');
    T.B(i).tEz=A';
    
    S=fgetl(fid);
    Ind=strfind(S,':');
    A=sscanf(S(Ind+1:end),'%e');
    T.B(i).nEx=A';

    S=fgetl(fid);
    Ind=strfind(S,':');
    A=sscanf(S(Ind+1:end),'%e');
    T.B(i).nEy=A';    

    S=fgetl(fid);
    Ind=strfind(S,':');
    A=sscanf(S(Ind+1:end),'%e');
    T.B(i).nEz=A'; 
    
    S=fgetl(fid);
    Ind=strfind(S,':');
    A=sscanf(S(Ind+1:end),'%e');
    T.B(i).sEx=A';

    S=fgetl(fid);
    Ind=strfind(S,':');
    A=sscanf(S(Ind+1:end),'%e');
    T.B(i).sEy=A';    

    S=fgetl(fid);
    Ind=strfind(S,':');
    A=sscanf(S(Ind+1:end),'%e');
    T.B(i).sEz=A';
    
    S=fgetl(fid);
    Ind=strfind(S,':');
    A=sscanf(S(Ind+1:end),'%e');
    T.B(i).ECtoR=A';
    
    S=fgetl(fid);
    Ind=strfind(S,':');
    A=sscanf(S(Ind+1:end),'%e');
    T.B(i).EAreaR=A';
    
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
    T.S(i).TtoC=A;
    
    S=fgetl(fid);
    Ind=strfind(S,':');
    A=sscanf(S(Ind+1:end),'%e');
    T.S(i).MCx=A';

    S=fgetl(fid);
    Ind=strfind(S,':');
    A=sscanf(S(Ind+1:end),'%e');
    T.S(i).MCy=A';
    
    S=fgetl(fid);
    Ind=strfind(S,':');
    A=sscanf(S(Ind+1:end),'%e');
    T.S(i).MCz=A';
    
    S=fgetl(fid);
    Ind=strfind(S,':');
    A=sscanf(S(Ind+1:end),'%e');
    T.S(i).CtoR=A';
    
    S=fgetl(fid);
    Ind=strfind(S,':');
    A=sscanf(S(Ind+1:end),'%e');
    T.S(i).PEx=A';

    S=fgetl(fid);
    Ind=strfind(S,':');
    A=sscanf(S(Ind+1:end),'%e');
    T.S(i).PEy=A';
    
    S=fgetl(fid);
    Ind=strfind(S,':');
    A=sscanf(S(Ind+1:end),'%e');
    T.S(i).PEz=A';

    S=fgetl(fid);
    Ind=strfind(S,':');
    A=sscanf(S(Ind+1:end),'%e');
    T.S(i).sEx=A';

    S=fgetl(fid);
    Ind=strfind(S,':');
    A=sscanf(S(Ind+1:end),'%e');
    T.S(i).sEy=A';
    
    S=fgetl(fid);
    Ind=strfind(S,':');
    A=sscanf(S(Ind+1:end),'%e');
    T.S(i).sEz=A';
    
    S=fgetl(fid);
    Ind=strfind(S,':');
    A=sscanf(S(Ind+1:end),'%e');
    T.S(i).ECtoR=A';
    
    S=fgetl(fid);
    Ind=strfind(S,':');
    A=sscanf(S(Ind+1:end),'%e');
    T.S(i).EAreaR=A';
    
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
