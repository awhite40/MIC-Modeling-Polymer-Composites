function [C_avg] = My_Abaqus_file_read( filename, n )
%UNTITLED3 Summary of this function goes here
%   n = number of elements on a side of a cube currently 21 
rawdata = importdata(filename,' ',1e8);
if size(rawdata,1) < 1000
   C_avg = NaN(3,3,3,3);
   Cs_avg = NaN(3,3,3,3);
   
else
%%
load('st1.mat')
%n=21;
line = 100;
E = cell((n^3)*8,9);
S = cell((n^3)*8,9);
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
E_star = zeros(size(E,1)/8,size(E,2)-1);
S1 = str2double(S);
S_star = zeros(size(S,1)/8,size(S,2)-1);
Emat = zeros(3,3,size(E,1)/8);
Smat = zeros(3,3,size(S,1)/8);
Savg = NaN(3,3);
Eavg = NaN(3,3);
C_avg = NaN(3,3,3,3);
Cs = NaN(3,3,3,3,n^3);
Cs_avg = NaN(3,3,3,3);
for ind =1:8:size(E,1)
    E_star(t,:) = ((E1(ind,2:9))+(E1(ind+1,2:9))+(E1(ind+2,2:9))+(E1(ind+3,2:9))...
        +(E1(ind+4,2:9))+(E1(ind+5,2:9))+(E1(ind+6,2:9))+(E1(ind+7,2:9)))/8;
    S_star(t,:) = ((S1(ind,2:9))+(S1(ind+1,2:9))+(S1(ind+2,2:9))+(S1(ind+3,2:9))...
        +(S1(ind+4,2:9))+(S1(ind+5,2:9))+(S1(ind+6,2:9))+(S1(ind+7,2:9)))/8;
    t=t+1;
end
Smat(1,1,:) = S_star(:,3);
Smat(1,2,:) = S_star(:,4);
Smat(2,1,:) = S_star(:,4);
Smat(1,3,:) = S_star(:,5);
Smat(3,1,:) = S_star(:,5);
Smat(2,2,:) = S_star(:,6);
Smat(2,3,:) = S_star(:,7);
Smat(3,2,:) = S_star(:,7);
Smat(3,3,:) = S_star(:,8);

        


Emat(1,1,:) = E_star(:,3);
Emat(1,2,:) = E_star(:,4);
Emat(2,1,:) = E_star(:,4);
Emat(1,3,:) = E_star(:,5);
Emat(3,1,:) = E_star(:,5);
Emat(2,2,:) = E_star(:,6);
Emat(2,3,:) = E_star(:,7);
Emat(3,2,:) = E_star(:,7);
Emat(3,3,:) = E_star(:,8);


i=1;
while i<=3
    j=1;
    while j<=3
        Savg(i,j) = mean(Smat(i,j,:));
        Eavg(i,j) = mean(Emat(i,j,:));
        if Eavg(i,j)<1e-10
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
        k=1;
        while k<=3
            l=1;
            while l<=3
                C_avg(i,j,k,l) = Savg(i,j)./Eavg(k,l);
                Cs(i,j,k,l,:) = Smat(i,k,:)./Emat(k,l,:);
                l=l+1;
            end
            k=k+1;
        end
        j=j+1;
    end
    i = i+1;
end

i=1;
while i<=3
    j=1;
    while j<=3
        k=1;
        while k<=3
            l=1;
            while l<=3
                Cs_avg(i,j,k,l) = mean(Cs(i,j,k,l,:));
                l=l+1;
            end
            k=k+1;
        end
        j=j+1;
    end
    i = i+1;
end
end
end