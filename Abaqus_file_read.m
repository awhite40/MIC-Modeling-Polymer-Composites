
[rawdata] = importdata('DiagFib4.dat',' ',1e8);
n = 21;
% st = rawdata(272,1);
% st2 = rawdata(74377,1);
load('st1.mat')
line = 1;
%E = cell((n^3)*8,9);
%S = cell((n^3)*8,9);
%%
while line <=size(rawdata,1)
    
    if strcmp(st,rawdata(line,1))
        chars = char(rawdata(line+3:line+2+(n^3)*8,1));
        for j = 1:size(chars,1)
            E(j,:) = strsplit(chars(j,:));
        end
        line = line +2+n^3*8;
    end
    
    if strcmp(st2,rawdata(line,1))
        chars = char(rawdata(line+3:line+2+(n^3)*8,1));
        for j = 1:size(chars,1)
            S(j,:) = strsplit(chars(j,:));
        end
        line = line +2+n^3*8;
    end
    
    line = line+1;
    
end
%%
E1 = str2double(E);
t=1;
E_star = zeros(size(E,1)/8,size(E,2)-3);
S1 = str2double(S);
S_star = zeros(size(S,1)/8,size(S,2)-3);
Emat = zeros(3,3,size(E,1)/8);
Smat = zeros(3,3,size(S,1)/8);
for ind =1:8:size(E,1)
    E_star(t,:) = ((E1(ind,4:9))+(E1(ind+1,4:9))+(E1(ind+2,4:9))+(E1(ind+3,4:9))...
        +(E1(ind+4,4:9))+(E1(ind+5,4:9))+(E1(ind+6,4:9))+(E1(ind+7,4:9)))/8;
    
    S_star(t,:) = ((S1(ind,4:9))+(S1(ind+1,4:9))+(S1(ind+2,4:9))+(S1(ind+3,4:9))...
        +(S1(ind+4,4:9))+(S1(ind+5,4:9))+(S1(ind+6,4:9))+(S1(ind+7,4:9)))/8;
    
    t=t+1;
end

Smat(1,1,:) = S_star(:,1);
Smat(1,2,:) = S_star(:,2);
Smat(2,1,:) = S_star(:,2);
Smat(1,3,:) = S_star(:,3);
Smat(3,1,:) = S_star(:,3);
Smat(2,2,:) = S_star(:,4);
Smat(2,3,:) = S_star(:,5);
Smat(3,2,:) = S_star(:,5);
Smat(3,3,:) = S_star(:,6);



Emat(1,1,:) = E_star(:,1);
Emat(1,2,:) = E_star(:,2);
Emat(2,1,:) = E_star(:,2);
Emat(1,3,:) = E_star(:,3);
Emat(3,1,:) = E_star(:,3);
Emat(2,2,:) = E_star(:,4);
Emat(2,3,:) = E_star(:,5);
Emat(3,2,:) = E_star(:,5);
Emat(3,3,:) = E_star(:,6);
%%
% Smat = ones(3,3,10);
% Emat = ones(3,3,10);
% 
% Smat(1,1,:)=10;
% Emat(1,2,1)= 5;


%%
i=1;
while i<=3
    j=1;
    while j<=3
        Savg(i,j) = mean(Smat(i,j,:));
        Eavg(i,j) = mean(Emat(i,j,:));
        if Eavg(i,j)<1e-8
            Eavg(i,j) =0;
        end
        j=j+1;
    end
    i=i+1;
end


i=1;
while i<=3
    j=1;
    while j<=3
        C_avg(i,j) = Savg(i,j)./Eavg(1,1);
            
        j=j+1;
    end
    i = i+1;
end
