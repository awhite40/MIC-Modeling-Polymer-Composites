clear all
close all
clc

%% Generate microstructure

s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s);

x=21; y=21; z=21;
H = 2; % Number of Phases
n = [60:5:120]/100; % Volume fraction
% n=1
filter = 5;

Iso1 = 3:2:15; % Dimension 1 (along x)
Iso2 = 1:2; % Dimension 2 (along y)
Iso3 = 1:2; % Dimension 3 (along z)
% Iso2 = 2;
% Iso3 = 2;

ct = 1;
for ii = 1:length(Iso1)
    for kk = 1:length(Iso2)
        for ll = 1:length(Iso3)
            for mm = 1:length(n)
                D = rand(x,y,z);
                [pp(:,:,:,ct), Vf(:,ct), f] = MsGenerator(D,H,n(mm),Iso1(ii),Iso2(kk),Iso3(ll),filter);
                ct = ct + 1;
%                 break
            end
%             break
        end
%         break
    end
%     break
end
% keyboard

ElongYFib = pp;
save ElongYFib.mat ElongYFib

%%

for ii =  1:size(pp,4)
    
    P = pp(:,:,:,ii);
    
    % Finding 1's and 0's
    dim = size(P);
    linInd = 1:prod(dim)';
    
    tmp = P(:);
    tmpEl = [linInd' tmp];
    %% Generate Input Files for random fibrous microstructrure
    
    % % % % % % % x=51; y=51; z=51;
    % % % % % % % D = rand(x,y,z);
    % % % % % % %
    % % % % % % % H = 2; % Number of Phases
    % % % % % % % n = 1.5; % Control volume fraction of phases
    % % % % % % % r = 0; % radial Grain size
    % % % % % % % Iso1 = 2; % Change in Grain size (stretch in 1-Dir)
    % % % % % % % Iso2 = 1;
    % % % % % % % P = MsGenerator(D,H,n,Iso1,Iso2,2,1);
    % % % % % % % if length(size(D)) < 3;
    % % % % % % %     pcolor(P);
    % % % % % % % else
    % % % % % % %     %     sliceomatic(P);
    % % % % % % % end
    % % % % % % %
    % % % % % % % % Finding 1's and 0's
    % % % % % % % dim = size(P);
    % % % % % % % linInd = 1:prod(dim)';
    % % % % % % %
    % % % % % % % tmp = P(:);
    % % % % % % % tmpEl = [linInd' tmp];
    % % % % % % %
    % % % % % % % % P(:) = reshape(tmpEl(:,2),dim);
    % % % % % % % keyboard
    
    %% Defining the dimension in abaqus
    
    nmfile = ['ElongYFib',num2str(ii) '.inp'];
    [SampXYZ] = ModelCoordinates(x,y,z,x,y,z);
    
    blockNo = 1;
    NodeNo1 = 1;
    ElementNo1 = 1;
    x0=0; y0=0; z0=0;
    Fine_int = 1;
    intx = x;
    inty = y;
    intz = z;
    Regions = 1;
    
    % Defining the most extreme nodes of the sample and its corresponding
    % coordinates
    
    setA=NodeNo1;
    setC=setA+(inty(1)*Fine_int);
    setB=setA+(inty(1)*Fine_int+1)*(intz(1)*Fine_int);
    setD=setA+(inty(1)*Fine_int+1)*(intz(1)*Fine_int)+(inty(1)*Fine_int);
    
    LayerNodes=(inty(1)*Fine_int+1)*(intz(1)*Fine_int)+(inty(1)*Fine_int)+1; % Not actual..this is with finer mesh
    
    setE=(intx(1)*Fine_int)*LayerNodes+setA;
    setG=setE+(inty(1)*Fine_int);
    setF=setE+(inty(1)*Fine_int+1)*(intz(1)*Fine_int);
    setH=setE+(inty(1)*Fine_int+1)*(intz(1)*Fine_int)+(inty(1)*Fine_int);
    
    SampXYZ=[setA SampXYZ(1,:,1);...
        setB SampXYZ(2,:,1);...
        setC SampXYZ(3,:,1);...
        setD SampXYZ(4,:,1);...
        setE SampXYZ(5,:,1);...
        setF SampXYZ(6,:,1);...
        setG SampXYZ(7,:,1);...
        setH SampXYZ(8,:,1)];
    
    %%
    temp=['*Heading'];
    dlmwrite(nmfile,temp,'delimiter','')
    temp=['** FE Simulation of Fibrous microstructure, By - Dipen Patel,  Date: ', date];
    
    dlmwrite(nmfile,temp,'delimiter','','-append')
    temp=['*Preprint, echo=NO, model=NO, history=NO, contact=NO'];
    dlmwrite(nmfile,temp,'delimiter','','-append')
    
    % 1] Print nodes on the corner with their coordinates
    temp = '***********************************************';
    dlmwrite(nmfile,temp,'delimiter','','-append')
    temp=['*node'];
    dlmwrite(nmfile,temp,'delimiter','','-append')
    dlmwrite(nmfile,SampXYZ,'delimiter',',','precision',12','-append')
    
    setAB_f=[intz(1) (inty(1)*Fine_int+1)*Fine_int];
    setCD_f=[intz(1) (inty(1)*Fine_int+1)*Fine_int];
    setABCD_f=[inty(1) Fine_int];
    
    setEF_f=[intz(1) (inty(1)*Fine_int+1)*Fine_int];
    setGH_f=[intz(1) (inty(1)*Fine_int+1)*Fine_int];
    setEFGH_f=[inty(1) Fine_int];
    
    setABCDEFGH_f=[intx(1) LayerNodes*Fine_int];
    
    % Define Element#1:
    Cube_1=setA;
    Cube_2=setA+Fine_int;
    Cube_3=setA+(inty(1)*Fine_int+2)*Fine_int;
    Cube_4=setA+(inty(1)*Fine_int+1)*Fine_int;
    Cube_5=(1*Fine_int)*LayerNodes+setA;
    Cube_6=(1*Fine_int)*LayerNodes+setA+Fine_int;
    Cube_7=Cube_5+(inty(1)*Fine_int+2)*Fine_int;
    Cube_8=Cube_5+(inty(1)*Fine_int+1)*Fine_int;
    
    Element1=[ElementNo1 Cube_1 Cube_2 Cube_3 Cube_4 Cube_5 Cube_6 Cube_7 Cube_8];
    
    Row_ELM_f=[inty(1) Fine_int Fine_int];
    Layer_ELM_f=[intz(1) setAB_f(2) (inty(1)*Fine_int)*Fine_int];
    Block_ELM_f=[intx(1) setABCDEFGH_f(2) (inty(1)*Fine_int)*(intz(1)*Fine_int)*Fine_int];
    
    ELGEN=[ElementNo1 Row_ELM_f Layer_ELM_f Block_ELM_f];
    
    [temp]=PrintFineMesh(...
        setA,setB,setC,setD,...
        setE,setF,setG,setH,...
        setAB_f, setCD_f,setEF_f,setGH_f,...
        setABCD_f,setEFGH_f,setABCDEFGH_f,...
        Element1,ELGEN,...
        blockNo,Regions,nmfile);
    
    %% Assigning material properties
    % tmpEl = [linInd' tmp];
    
    tmp  = unique(tmpEl(:,2));
    
    mat1 = ismember(tmpEl(:,2),tmp(1,1));
    mat2 = ismember(tmpEl(:,2),tmp(2,1));
    
    mat1inx = find(mat1 ==1);
    mat2inx = find(mat2 ==1);
    
    
    if mod(length(mat1inx),16) == 0
        mat1inx = reshape(mat1inx,16,[]);
        temp=['*Elset, elset=Matx'];
        dlmwrite(nmfile,temp,'delimiter','','-append')
        dlmwrite(nmfile,mat1inx(1:end,:),'precision',12,'delimiter',',','-append')
        %     dlmwrite(nmfile,mat1inx(end,1:rim),'precision',12,'delimiter',',','-append')
    else
        Orim = mod(length(mat1inx),16);
        rim = 16-Orim;
        mat1inx(end+rim) = 0;
        mat1inx = reshape(mat1inx,16,[])';
        
        temp=['*Elset, elset=Matx'];
        dlmwrite(nmfile,temp,'delimiter','','-append')
        dlmwrite(nmfile,mat1inx(1:end-1,:),'precision',12,'delimiter',',','-append')
        dlmwrite(nmfile,mat1inx(end,1:Orim),'precision',12,'delimiter',',','-append')
    end
    
    if mod(length(mat2inx),16) == 0
        mat2inx = reshape(mat2inx,16,[]);
        temp=['*Elset, elset=Fib'];
        dlmwrite(nmfile,temp,'delimiter','','-append')
        dlmwrite(nmfile,mat2inx(1:end,:),'precision',12,'delimiter',',','-append')
        %     dlmwrite(nmfile,mat2inx(end,1:rim),'precision',12,'delimiter',',','-append')
    else
        Orim = mod(length(mat2inx),16);
        rim = 16-Orim;
        mat2inx(end+rim) = 0;
        mat2inx = reshape(mat2inx,16,[])';
        temp=['*Elset, elset=Fib'];
        dlmwrite(nmfile,temp,'delimiter','','-append')
        dlmwrite(nmfile,mat2inx(1:end-1,:),'precision',12,'delimiter',',','-append')
        dlmwrite(nmfile,mat2inx(end,1:Orim),'precision',12,'delimiter',',','-append')
    end
    
    %% Preparing the nodes for periodic boundary conditions
    % n1minus
    temp = (((setA + (x+1)) + 1 ):1:(setABCD_f(:,1)*2)+1 )';
    for ii = 1:x-2
        tmp1 = temp(:,ii) + (y+1);
        temp(:,ii+1) = tmp1;
    end
    temp = temp(:);
    
    if mod(length(temp),15) == 0
        temp = reshape(temp,15,[]);
        tmp=['*Nset, Nset=n1minus'];
        dlmwrite(nmfile,tmp,'delimiter','','-append')
        dlmwrite(nmfile,temp(1:end,:),'precision',12,'delimiter',',','-append')
        %     dlmwrite(nmfile,temp(end,1:rim),'precision',12,'delimiter',',','-append')
    else
        Orim = mod(length(temp),15);
        rim = 15-Orim;
        temp(end+rim) = 0;
        temp = reshape(temp,15,[])';
        tmp=['*Nset, Nset=n1minus'];
        dlmwrite(nmfile,tmp,'delimiter','','-append')
        dlmwrite(nmfile,temp(1:end-1,:),'precision',12,'delimiter',',','-append')
        dlmwrite(nmfile,temp(end,1:Orim),'precision',12,'delimiter',',','-append')
    end
    
    % n1plus
    temp = (((setE + (x+1)) + 1 ):1:(setABCD_f(:,1)*2)+setE )';
    for ii = 1:x-2
        tmp1 = temp(:,ii) + (y+1);
        temp(:,ii+1) = tmp1;
    end
    temp = temp(:);
    
    if mod(length(temp),16) == 0
        temp = reshape(temp,16,[])';
        tmp=['*Nset, Nset=n1plus'];
        dlmwrite(nmfile,tmp,'delimiter','','-append')
        dlmwrite(nmfile,temp(1:end,:),'precision',12,'delimiter',',','-append')
        %     dlmwrite(nmfile,temp(end,1:rim),'precision',12,'delimiter',',','-append')
    else
        Orim = mod(length(temp),16);
        rim = 16-Orim;
        temp(end+rim) = 0;
        temp = reshape(temp,16,[])';
        tmp=['*Nset, Nset=n1plus'];
        dlmwrite(nmfile,tmp,'delimiter','','-append')
        dlmwrite(nmfile,temp(1:end-1,:),'precision',12,'delimiter',',','-append')
        dlmwrite(nmfile,temp(end,1:Orim),'precision',12,'delimiter',',','-append')
    end
    
    %% n2minus
    temp = (((setD+1)+setAB_f(:,1)+1):setAB_f(:,2):((setD+1)+setAB_f(:,1)+1)+(setAB_f(:,2)*(y-2)))';
    for ii = 1:z-2
        tmp1 = temp(:,ii) + LayerNodes;
        temp(:,ii+1) = tmp1;
    end
    temp = temp(:);
    if mod(length(temp),16) == 0
        temp = reshape(temp,16,[])';
        tmp=['*Nset, Nset=n2minus'];
        dlmwrite(nmfile,tmp,'delimiter','','-append')
        dlmwrite(nmfile,temp(1:end,:),'precision',12,'delimiter',',','-append')
        %     dlmwrite(nmfile,temp(end,1:rim),'precision',12,'delimiter',',','-append')
    else
        Orim = mod(length(temp),16);
        rim = 16-Orim;
        temp(end+rim) = 0;
        temp = reshape(temp,16,[])';
        tmp=['*Nset, Nset=n2minus'];
        dlmwrite(nmfile,tmp,'delimiter','','-append')
        dlmwrite(nmfile,temp(1:end-1,:),'precision',12,'delimiter',',','-append')
        dlmwrite(nmfile,temp(end,1:Orim),'precision',12,'delimiter',',','-append')
    end
    
    % n2plus
    temp = ((((setD+1)+setAB_f(:,1)+1)+(y-1)+1):setAB_f(:,2):(setD+1)+setAB_f(:,1)+(setAB_f(:,2)*(z-1)))';
    for ii = 1:z-2
        tmp1 = temp(:,ii) + LayerNodes;
        temp(:,ii+1) = tmp1;
    end
    temp = temp(:);
    if mod(length(temp),16) == 0
        temp = reshape(temp,16,[])';
        tmp=['*Nset, Nset=n2plus'];
        dlmwrite(nmfile,tmp,'delimiter','','-append')
        dlmwrite(nmfile,temp(1:end,:),'precision',12,'delimiter',',','-append')
        %     dlmwrite(nmfile,temp(end,1:rim),'precision',12,'delimiter',',','-append')
    else
        Orim = mod(length(temp),16);
        rim = 16-Orim;
        temp(end+rim) = 0;
        temp = reshape(temp,16,[])';
        tmp=['*Nset, Nset=n2plus'];
        dlmwrite(nmfile,tmp,'delimiter','','-append')
        dlmwrite(nmfile,temp(1:end-1,:),'precision',12,'delimiter',',','-append')
        dlmwrite(nmfile,temp(end,1:Orim),'precision',12,'delimiter',',','-append')
    end
    %% n3minus
    temp = ((setD+2):1:(setD)+(inty(1)*Fine_int))';
    for ii = 1:x-2
        tmp1 = temp(:,ii) + LayerNodes;
        temp(:,ii+1) = tmp1;
    end
    temp = temp(:);
    
    
    if mod(length(temp),16) == 0
        temp = reshape(temp,16,[])';
        tmp=['*Nset, Nset=n3minus'];
        dlmwrite(nmfile,tmp,'delimiter','','-append')
        dlmwrite(nmfile,temp(1:end,:),'precision',12,'delimiter',',','-append')
        %     dlmwrite(nmfile,temp(end,1:rim),'precision',12,'delimiter',',','-append')
    else
        Orim = mod(length(temp),16);
        rim = 16-Orim;
        temp(end+rim) = 0;
        temp = reshape(temp,16,[])';
        tmp=['*Nset, Nset=n3minus'];
        dlmwrite(nmfile,tmp,'delimiter','','-append')
        dlmwrite(nmfile,temp(1:end-1,:),'precision',12,'delimiter',',','-append')
        dlmwrite(nmfile,temp(end,1:Orim),'precision',12,'delimiter',',','-append')
    end
    
    %n3plus
    temp = ((setD+2) + (z*setAB_f(:,2)):1:((setD+2) + (z*setAB_f(:,2)))+((x-2)*setABCD_f(:,2)))';
    for ii = 1:x-2
        tmp1 = temp(:,ii) + LayerNodes;
        temp(:,ii+1) = tmp1;
    end
    temp = temp(:);
    
    if mod(length(temp),16) == 0
        temp = reshape(temp,16,[])';
        tmp=['*Nset, Nset=n3plus'];
        dlmwrite(nmfile,tmp,'delimiter','','-append')
        dlmwrite(nmfile,temp(1:end,:),'precision',12,'delimiter',',','-append')
        %     dlmwrite(nmfile,temp(end,1:rim),'precision',12,'delimiter',',','-append')
    else
        Orim = mod(length(temp),16);
        rim = 16-Orim;
        temp(end+rim) = 0;
        temp = reshape(temp,16,[])';
        tmp=['*Nset, Nset=n3plus'];
        dlmwrite(nmfile,tmp,'delimiter','','-append')
        dlmwrite(nmfile,temp(1:end-1,:),'precision',12,'delimiter',',','-append')
        dlmwrite(nmfile,temp(end,1:Orim),'precision',12,'delimiter',',','-append')
    end
    
    %% Preparing Edges
    
    temp = (setA+1:setC-1)';
    if mod(length(temp),16) == 0
        temp = reshape(temp,16,[])';
        tmp=['*NSET, NSET=n3minus_n1minus'];
        dlmwrite(nmfile,tmp,'delimiter','','-append')
        dlmwrite(nmfile,temp(1:end,:),'precision',12,'delimiter',',','-append')
        %     dlmwrite(nmfile,temp(end,1:rim),'precision',12,'delimiter',',','-append')
    else
        Orim = mod(length(temp),16);
        rim = 16-Orim;
        temp(end+rim) = 0;
        temp = reshape(temp,16,[])';
        tmp=['*NSET, NSET=n3minus_n1minus'];
        dlmwrite(nmfile,tmp,'delimiter','','-append')
        dlmwrite(nmfile,temp(1:end-1,:),'precision',12,'delimiter',',','-append')
        dlmwrite(nmfile,temp(end,1:Orim),'precision',12,'delimiter',',','-append')
    end
    
    temp = (setB+1:setD-1)';
    if mod(length(temp),16) == 0
        temp = reshape(temp,16,[])';
        tmp=['*NSET, NSET=n3plus_n1minus'];
        dlmwrite(nmfile,tmp,'delimiter','','-append')
        dlmwrite(nmfile,temp(1:end,:),'precision',12,'delimiter',',','-append')
        %     dlmwrite(nmfile,temp(end,1:rim),'precision',12,'delimiter',',','-append')
    else
        Orim = mod(length(temp),16);
        rim = 16-Orim;
        temp(end+rim) = 0;
        temp = reshape(temp,16,[])';
        tmp=['*NSET, NSET=n3plus_n1minus'];
        dlmwrite(nmfile,tmp,'delimiter','','-append')
        dlmwrite(nmfile,temp(1:end-1,:),'precision',12,'delimiter',',','-append')
        dlmwrite(nmfile,temp(end,1:Orim),'precision',12,'delimiter',',','-append')
    end
    
    temp = (setE+1:setG-1)';
    if mod(length(temp),16) == 0
        temp = reshape(temp,16,[])';
        tmp=['*NSET, NSET=n3minus_n1plus'];
        dlmwrite(nmfile,tmp,'delimiter','','-append')
        dlmwrite(nmfile,temp(1:end,:),'precision',12,'delimiter',',','-append')
        %     dlmwrite(nmfile,temp(end,1:rim),'precision',12,'delimiter',',','-append')
    else
        Orim = mod(length(temp),16);
        rim = 16-Orim;
        temp(end+rim) = 0;
        temp = reshape(temp,16,[])';
        tmp=['*NSET, NSET=n3minus_n1plus'];
        dlmwrite(nmfile,tmp,'delimiter','','-append')
        dlmwrite(nmfile,temp(1:end-1,:),'precision',12,'delimiter',',','-append')
        dlmwrite(nmfile,temp(end,1:Orim),'precision',12,'delimiter',',','-append')
    end
    
    
    temp = (setF+1:setH-1)';
    if mod(length(temp),16) == 0
        temp = reshape(temp,16,[])';
        tmp=['*NSET, NSET=n3plus_n1plus'];
        dlmwrite(nmfile,tmp,'delimiter','','-append')
        dlmwrite(nmfile,temp(1:end,:),'precision',12,'delimiter',',','-append')
        %     dlmwrite(nmfile,temp(end,1:rim),'precision',12,'delimiter',',','-append')
    else
        Orim = mod(length(temp),16);
        rim = 16-Orim;
        temp(end+rim) = 0;
        temp = reshape(temp,16,[])';
        tmp=['*NSET, NSET=n3plus_n1plus'];
        dlmwrite(nmfile,tmp,'delimiter','','-append')
        dlmwrite(nmfile,temp(1:end-1,:),'precision',12,'delimiter',',','-append')
        dlmwrite(nmfile,temp(end,1:Orim),'precision',12,'delimiter',',','-append')
    end
    
    temp = (setA+setAB_f(:,1)+1:setAB_f(:,2):setB-1)';
    if mod(length(temp),16) == 0
        temp = reshape(temp,16,[])';
        tmp=['*NSET, NSET=n1minus_n2minus'];
        dlmwrite(nmfile,tmp,'delimiter','','-append')
        dlmwrite(nmfile,temp(1:end,:),'precision',12,'delimiter',',','-append')
        %     dlmwrite(nmfile,temp(end,1:rim),'precision',12,'delimiter',',','-append')
    else
        Orim = mod(length(temp),16);
        rim = 16-Orim;
        temp(end+rim) = 0;
        temp = reshape(temp,16,[])';
        tmp=['*NSET, NSET=n1minus_n2minus'];
        dlmwrite(nmfile,tmp,'delimiter','','-append')
        dlmwrite(nmfile,temp(1:end-1,:),'precision',12,'delimiter',',','-append')
        dlmwrite(nmfile,temp(end,1:Orim),'precision',12,'delimiter',',','-append')
    end
    
    
    temp = (setC+setAB_f(:,1)+1:setAB_f(:,2):setD-1)';
    if mod(length(temp),16) == 0
        temp = reshape(temp,16,[])';
        tmp=['*NSET, NSET=n1minus_n2plus'];
        dlmwrite(nmfile,tmp,'delimiter','','-append')
        dlmwrite(nmfile,temp(1:end,:),'precision',12,'delimiter',',','-append')
        %     dlmwrite(nmfile,temp(end,1:rim),'precision',12,'delimiter',',','-append')
    else
        Orim = mod(length(temp),16);
        rim = 16-Orim;
        temp(end+rim) = 0;
        temp = reshape(temp,16,[])';
        tmp=['*NSET, NSET=n1minus_n2plus'];
        dlmwrite(nmfile,tmp,'delimiter','','-append')
        dlmwrite(nmfile,temp(1:end-1,:),'precision',12,'delimiter',',','-append')
        dlmwrite(nmfile,temp(end,1:Orim),'precision',12,'delimiter',',','-append')
    end
    
    temp = (setE+setAB_f(:,1)+1:setAB_f(:,2):setF-1)';
    if mod(length(temp),16) == 0
        temp = reshape(temp,16,[])';
        tmp=['*NSET, NSET=n1plus_n2minus'];
        dlmwrite(nmfile,tmp,'delimiter','','-append')
        dlmwrite(nmfile,temp(1:end,:),'precision',12,'delimiter',',','-append')
        %     dlmwrite(nmfile,temp(end,1:rim),'precision',12,'delimiter',',','-append')
    else
        Orim = mod(length(temp),16);
        rim = 16-Orim;
        temp(end+rim) = 0;
        temp = reshape(temp,16,[])';
        tmp=['*NSET, NSET=n1plus_n2minus'];
        dlmwrite(nmfile,tmp,'delimiter','','-append')
        dlmwrite(nmfile,temp(1:end-1,:),'precision',12,'delimiter',',','-append')
        dlmwrite(nmfile,temp(end,1:Orim),'precision',12,'delimiter',',','-append')
    end
    
    temp =(setG+setAB_f(:,1)+1:setAB_f(:,2):setH-1)';
    if mod(length(temp),16) == 0
        temp = reshape(temp,16,[])';
        tmp=['*NSET, NSET=n1plus_n2plus'];
        dlmwrite(nmfile,tmp,'delimiter','','-append')
        dlmwrite(nmfile,temp(1:end,:),'precision',12,'delimiter',',','-append')
        %     dlmwrite(nmfile,temp(end,1:rim),'precision',12,'delimiter',',','-append')
    else
        Orim = mod(length(temp),16);
        rim = 16-Orim;
        temp(end+rim) = 0;
        temp = reshape(temp,16,[])';
        tmp=['*NSET, NSET=n1plus_n2plus'];
        dlmwrite(nmfile,tmp,'delimiter','','-append')
        dlmwrite(nmfile,temp(1:end-1,:),'precision',12,'delimiter',',','-append')
        dlmwrite(nmfile,temp(end,1:Orim),'precision',12,'delimiter',',','-append')
    end
    
    
    temp= (setA+LayerNodes:LayerNodes:setA+((x-1)*LayerNodes))';
    if mod(length(temp),16) == 0
        temp = reshape(temp,16,[])';
        tmp=['*NSET, NSET=n2minus_n3minus'];
        dlmwrite(nmfile,tmp,'delimiter','','-append')
        dlmwrite(nmfile,temp(1:end,:),'precision',12,'delimiter',',','-append')
        %     dlmwrite(nmfile,temp(end,1:rim),'precision',12,'delimiter',',','-append')
    else
        Orim = mod(length(temp),16);
        rim = 16-Orim;
        temp(end+rim) = 0;
        temp = reshape(temp,16,[])';
        tmp=['*NSET, NSET=n2minus_n3minus'];
        dlmwrite(nmfile,tmp,'delimiter','','-append')
        dlmwrite(nmfile,temp(1:end-1,:),'precision',12,'delimiter',',','-append')
        dlmwrite(nmfile,temp(end,1:Orim),'precision',12,'delimiter',',','-append')
    end
    
    temp = (setC+LayerNodes:LayerNodes:setC+((x-1)*LayerNodes))';
    if mod(length(temp),16) == 0
        temp = reshape(temp,16,[])';
        tmp=['*NSET, NSET=n2plus_n3minus'];
        dlmwrite(nmfile,tmp,'delimiter','','-append')
        dlmwrite(nmfile,temp(1:end,:),'precision',12,'delimiter',',','-append')
        %     dlmwrite(nmfile,temp(end,1:rim),'precision',12,'delimiter',',','-append')
    else
        Orim = mod(length(temp),16);
        rim = 16-Orim;
        temp(end+rim) = 0;
        temp = reshape(temp,16,[])';
        tmp=['*NSET, NSET=n2plus_n3minus'];
        dlmwrite(nmfile,tmp,'delimiter','','-append')
        dlmwrite(nmfile,temp(1:end-1,:),'precision',12,'delimiter',',','-append')
        dlmwrite(nmfile,temp(end,1:Orim),'precision',12,'delimiter',',','-append')
    end
    
    temp = (setB+LayerNodes:LayerNodes:setB+((x-1)*LayerNodes))';
    if mod(length(temp),16) == 0
        temp = reshape(temp,16,[])';
        tmp=['*NSET, NSET=n2minus_n3plus'];
        dlmwrite(nmfile,tmp,'delimiter','','-append')
        dlmwrite(nmfile,temp(1:end,:),'precision',12,'delimiter',',','-append')
        %     dlmwrite(nmfile,temp(end,1:rim),'precision',12,'delimiter',',','-append')
    else
        Orim = mod(length(temp),16);
        rim = 16-Orim;
        temp(end+rim) = 0;
        temp = reshape(temp,16,[])';
        tmp=['*NSET, NSET=n2minus_n3plus'];
        dlmwrite(nmfile,tmp,'delimiter','','-append')
        dlmwrite(nmfile,temp(1:end-1,:),'precision',12,'delimiter',',','-append')
        dlmwrite(nmfile,temp(end,1:Orim),'precision',12,'delimiter',',','-append')
    end
    
    temp = setD+LayerNodes:LayerNodes:setD+((x-1)*LayerNodes);
    if mod(length(temp),16) == 0
        temp = reshape(temp,16,[])';
        tmp=['*NSET, NSET=n2plus_n3plus'];
        dlmwrite(nmfile,tmp,'delimiter','','-append')
        dlmwrite(nmfile,temp(1:end,:),'precision',12,'delimiter',',','-append')
        %     dlmwrite(nmfile,temp(end,1:rim),'precision',12,'delimiter',',','-append')
    else
        Orim = mod(length(temp),16);
        rim = 16-Orim;
        temp(end+rim) = 0;
        temp = reshape(temp,16,[])';
        tmp=['*NSET, NSET=n2plus_n3plus'];
        dlmwrite(nmfile,tmp,'delimiter','','-append')
        dlmwrite(nmfile,temp(1:end-1,:),'precision',12,'delimiter',',','-append')
        dlmwrite(nmfile,temp(end,1:Orim),'precision',12,'delimiter',',','-append')
    end
    
    %% Applying periodic boundary conditions and displacement on one of the faces.
    
    tmp='*Equation';
    dlmwrite(nmfile,tmp,'delimiter','','-append');
    tmp='3';
    dlmwrite(nmfile,tmp,'delimiter','','-append')
    tmp=['n1plus, 1, 1, n1minus, 1, -1, ' num2str(setH) ', 1, -1'];
    dlmwrite(nmfile,tmp,'precision',12,'delimiter','','-append')
    
    tmp='3';
    dlmwrite(nmfile,tmp,'delimiter','','-append');
    tmp=['n1plus_n2plus, 1, 1, n1minus_n2plus, 1, -1,' num2str(setH) ', 1, -1'];
    dlmwrite(nmfile,tmp,'precision',12,'delimiter','','-append')
    
    tmp='3';
    dlmwrite(nmfile,tmp,'delimiter','','-append');
    tmp=['n1plus_n2minus, 1, 1, n1minus_n2minus, 1, -1,' num2str(setH) ', 1, -1'];
    dlmwrite(nmfile,tmp,'precision',12,'delimiter','','-append')
    
    tmp='3';
    dlmwrite(nmfile,tmp,'delimiter','','-append');
    tmp=['n3plus_n1plus, 1, 1, n3plus_n1minus, 1, -1,' num2str(setH) ', 1, -1'];
    dlmwrite(nmfile,tmp,'precision',12,'delimiter','','-append')
    
    tmp='3';
    dlmwrite(nmfile,tmp,'delimiter','','-append');
    tmp=['n3minus_n1plus, 1, 1, n3minus_n1minus, 1, -1,' num2str(setH) ', 1, -1'];
    dlmwrite(nmfile,tmp,'precision',12,'delimiter','','-append')
    
    tmp=['*Include, Input=condPBC.txt'];
    dlmwrite(nmfile,tmp,'delimiter','','-append')
    
    %% Applying Material properties
    
    tmp=['*Include, Input=MatProp.txt'];
    dlmwrite(nmfile,tmp,'delimiter','','-append')
    
    %% Applying displacement to one of the faces
    
    tmp='*Step, name=Step-1';
    dlmwrite(nmfile,tmp,'delimiter','','-append')
    tmp='*Static';
    dlmwrite(nmfile,tmp,'delimiter','','-append')
    tmp=[1., 1., 1e-05, 1.];
    dlmwrite(nmfile,tmp,'precision',12,'delimiter',',','-append')
    
    tmp='*Boundary';
    dlmwrite(nmfile,tmp,'delimiter','','-append')
    tmp=[setA,1,3,0];
    dlmwrite(nmfile,tmp,'precision',12,'delimiter',',','-append')
    tmp=[setB,1,3,0];
    dlmwrite(nmfile,tmp,'precision',12,'delimiter',',','-append')
    tmp=[setC,1,3,0];
    dlmwrite(nmfile,tmp,'precision',12,'delimiter',',','-append')
    tmp=[setD,1,3,0];
    dlmwrite(nmfile,tmp,'precision',12,'delimiter',',','-append')
    tmp=[setE,1,1,0.02];
    dlmwrite(nmfile,tmp,'precision',12,'delimiter',',','-append')
    tmp=[setF,1,1,0.02];
    dlmwrite(nmfile,tmp,'precision',12,'delimiter',',','-append')
    tmp=[setG,1,1,0.02];
    dlmwrite(nmfile,tmp,'precision',12,'delimiter',',','-append')
    tmp=[setH,1,1,0.02];
    dlmwrite(nmfile,tmp,'precision',12,'delimiter',',','-append')
    tmp=[setE,2,3,0];
    dlmwrite(nmfile,tmp,'precision',12,'delimiter',',','-append')
    tmp=[setF,2,3,0];
    dlmwrite(nmfile,tmp,'precision',12,'delimiter',',','-append')
    tmp=[setG,2,3,0];
    dlmwrite(nmfile,tmp,'precision',12,'delimiter',',','-append')
    tmp=[setH,2,3,0];
    dlmwrite(nmfile,tmp,'precision',12,'delimiter',',','-append')
    
    %% Prescribing the output (stress,strain)
    tmp=['*Include, Input=resOut.txt'];
    dlmwrite(nmfile,tmp,'delimiter','','-append')
    
    % % % tmp=['** OUTPUT REQUESTS'];
    % % % dlmwrite(nmfile,tmp,'delimiter','','-append')
    % % % tmp=['**'];
    % % % dlmwrite(nmfile,tmp,'delimiter','','-append')
    % % % tmp=['*output, field, frequency=0'];
    % % % dlmwrite(nmfile,tmp,'delimiter','','-append')
    % % % tmp=['**'];
    % % % dlmwrite(nmfile,tmp,'delimiter','','-append')
    % % % tmp=['*output, history, frequency=0'];
    % % % dlmwrite(nmfile,tmp,'delimiter','','-append')
    % % % tmp=['** '];
    % % % dlmwrite(nmfile,tmp,'delimiter','','-append')
    % % % tmp=['*el print, summary=no, totals=yes'];
    % % % dlmwrite(nmfile,tmp,'delimiter','','-append')
    % % % tmp=['E, S'];
    % % % dlmwrite(nmfile,tmp,'delimiter','','-append')
    % % % tmp=['**'];
    % % % dlmwrite(nmfile,tmp,'delimiter','','-append')
    % % % tmp=['*End Step'];
    % % % dlmwrite(nmfile,tmp,'delimiter','','-append')
    
end
