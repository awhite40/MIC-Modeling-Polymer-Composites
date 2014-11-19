[rawdata] = importdata('tst.dat',' ',1e8);
%FID = fopen('tst.dat');
%formatspec = '%f %f %f %f %f %f';
%fm2 = '%s %s %s %s %s %s %s %s';
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
chars = str2num(E(1,:))
E_star(1,:) = mean(char(E(1:8,:)));

