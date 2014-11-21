
[rawdata] = importdata('tst.dat',' ',1e8);
n = 21;
st = rawdata(268,1);
st2 = rawdata(74373,1);
line = 1;
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
E_star = zeros(size(E,1)/8,size(E,2)-3);
S1 = str2double(S);
S_star = zeros(size(S,1)/8,size(S,2)-1);
for ind =1:8:size(E,1)
    E_star(t,:) = ((E1(ind,4:9))+(E1(ind+1,4:9))+(E1(ind+2,4:9))+(E1(ind+3,4:9))...
        +(E1(ind+4,4:9))+(E1(ind+5,4:9))+(E1(ind+6,4:9))+(E1(ind+7,4:9)))/8;
    S_star(t,:) = ((S1(ind,2:9))+(S1(ind+1,2:9))+(S1(ind+2,2:9))+(S1(ind+3,2:9))...
        +(S1(ind+4,2:9))+(S1(ind+5,2:9))+(S1(ind+6,2:9))+(S1(ind+7,2:9)))/8;
    t=t+1;
end
