% 3D Modified Nanoindentation Model
% Author:- Dipen Patel
% Version # 1
% Please email me if you find any error


function [temp]=PrintFineMesh(...
    setA,setB,setC,setD,...
    setE,setF,setG,setH,...
    setAB, setCD,setEF,setGH,...
    setABCD,setEFGH,setABCDEFGH,...
    Element1,ELGEN,...
    blockNo,Regions,Name_input)

fil_nm = [Name_input]

temp=['***********************************************'];
dlmwrite(fil_nm,temp,'delimiter','','-append')

temp=['*NSET, NSET=A',num2str(blockNo)];
dlmwrite(fil_nm,temp,'delimiter','','-append')
dlmwrite(fil_nm,setA,'precision',12,'delimiter','','-append')


temp=['*NSET, NSET=B',num2str(blockNo)];
dlmwrite(fil_nm,temp,'delimiter','','-append')
dlmwrite(fil_nm,setB,'precision',12,'delimiter','','-append')


temp=['*NSET, NSET=C',num2str(blockNo)];
dlmwrite(fil_nm,temp,'delimiter','','-append')
dlmwrite(fil_nm,setC,'precision',12,'delimiter','','-append')


temp=['*NSET, NSET=D',num2str(blockNo)];
dlmwrite(fil_nm,temp,'delimiter','','-append')
dlmwrite(fil_nm,setD,'precision',12,'delimiter','','-append')


temp=['*NSET, NSET=E',num2str(blockNo)];
dlmwrite(fil_nm,temp,'delimiter','','-append')
dlmwrite(fil_nm,setE,'precision',12,'delimiter','','-append')


temp=['*NSET, NSET=F',num2str(blockNo)];
dlmwrite(fil_nm,temp,'delimiter','','-append')
dlmwrite(fil_nm,setF,'precision',12,'delimiter','','-append')


temp=['*NSET, NSET=G',num2str(blockNo)];
dlmwrite(fil_nm,temp,'delimiter','','-append')
dlmwrite(fil_nm,setG,'precision',12,'delimiter','','-append')


temp=['*NSET, NSET=H',num2str(blockNo)];
dlmwrite(fil_nm,temp,'delimiter','','-append')
dlmwrite(fil_nm,setH,'precision',12,'delimiter','','-append')


%% 2] Create nodes b/w edges
temp=['***********************************************'];
dlmwrite(fil_nm,temp,'delimiter','','-append')

temp=['*Nfill, Nset=AB',num2str(blockNo)];
dlmwrite(fil_nm,temp,'delimiter','','-append')
temp=['A' num2str(blockNo) ', B' num2str(blockNo) ', ' num2str(setAB(1)) ', ' num2str(setAB(2))];
dlmwrite(fil_nm,temp,'precision',12,'delimiter','','-append')

temp=['*Nfill, Nset=CD'  num2str(blockNo)];
dlmwrite(fil_nm,temp,'delimiter','','-append')
temp=['C' num2str(blockNo) ', D' num2str(blockNo) ', ' num2str(setCD(1)) ', ' num2str(setCD(2))];
dlmwrite(fil_nm,temp,'precision',12,'delimiter','','-append')


temp=['*Nfill, Nset=ABCD'  num2str(blockNo)];
dlmwrite(fil_nm,temp,'delimiter','','-append')
temp=['AB' num2str(blockNo) ', CD' num2str(blockNo) ', ' num2str(setABCD(1)) ', ' num2str(setABCD(2))];
dlmwrite(fil_nm,temp,'precision',12,'delimiter','','-append')

temp=['*Nfill, Nset=EF'  num2str(blockNo)];
dlmwrite(fil_nm,temp,'delimiter','','-append')
temp=['E' num2str(blockNo) ', F' num2str(blockNo) ', ' num2str(setEF(1)) ', ' num2str(setEF(2))];
dlmwrite(fil_nm,temp,'precision',12,'delimiter','','-append')

temp=['*Nfill, Nset=GH'  num2str(blockNo)];
dlmwrite(fil_nm,temp,'delimiter','','-append')
temp=['G' num2str(blockNo) ', H' num2str(blockNo) ', ' num2str(setGH(1)) ', ' num2str(setGH(2))];
dlmwrite(fil_nm,temp,'precision',12,'delimiter','','-append')

temp=['*Nfill, Nset=EFGH'  num2str(blockNo)];
dlmwrite(fil_nm,temp,'delimiter','','-append')
temp=['EF' num2str(blockNo) ', GH' num2str(blockNo) ', ' num2str(setEFGH(1)) ', ' num2str(setEFGH(2))];
dlmwrite(fil_nm,temp,'precision',12,'delimiter','','-append')

temp=['*Nfill, Nset=ABCDEFGH'  num2str(blockNo)];
dlmwrite(fil_nm,temp,'delimiter','','-append')
temp=['ABCD' num2str(blockNo) ', EFGH' num2str(blockNo) ', ' num2str(setABCDEFGH(1)) ', ' num2str(setABCDEFGH(2))];
dlmwrite(fil_nm,temp,'precision',12,'delimiter','','-append')


%%  Define the top surface of the finer mesh using Node-based approach 
%   See Fun_Nanoind_SampleMesh.m for Element-based approach
% if blockNo==Regions
% 
%     temp=['***********************************************'];
%     dlmwrite(fil_nm,temp,'delimiter','','-append')
%     
%     temp=['*Nfill, NSET=BD',num2str(blockNo)];
%     dlmwrite(fil_nm,temp,'delimiter','','-append')
%     temp=['B' num2str(blockNo) ', D' num2str(blockNo) ', ' num2str(setABCD(1)) ', ' num2str(setABCD(2))];
%     dlmwrite(fil_nm,temp,'precision',12,'delimiter','','-append')
%     
%     temp=['*Nfill, NSET=FH',num2str(blockNo)];
%     dlmwrite(fil_nm,temp,'delimiter','','-append')
%     temp=['F' num2str(blockNo) ', H' num2str(blockNo) ', ' num2str(setABCD(1)) ', ' num2str(setABCD(2))];
%     dlmwrite(fil_nm,temp,'precision',12,'delimiter','','-append')
%     
%     
%     temp=['*Nfill, NSET=BDFH',num2str(blockNo)];
%     dlmwrite(fil_nm,temp,'delimiter','','-append')
%     temp=['BD' num2str(blockNo) ', FH' num2str(blockNo) ', ' num2str(setABCDEFGH(1)) ', ' num2str(setABCDEFGH(2))];
%     dlmwrite(fil_nm,temp,'precision',12,'delimiter','','-append')
%     
%     
%     temp=['***********************************************'];
%     dlmwrite(fil_nm,temp,'delimiter','','-append')
%     temp=['*Surface, Name=TopfineSurf_NODEbased, Type=Node'];
%     dlmwrite(fil_nm,temp,'delimiter','','-append')
%     temp=['BDFH',num2str(blockNo)];
%     dlmwrite(fil_nm,temp,'delimiter','','-append')
%        
% end

%%  Define the Bottom surface of the sample using Node-based approach:
% if blockNo==1
%     temp=['***********************************************'];
%     dlmwrite(fil_nm,temp,'delimiter','','-append')
%     temp=['** Bottom Face Node Set'];
%     dlmwrite(fil_nm,temp,'delimiter','','-append')
%     
%     temp=['*Nfill, NSET=AC',num2str(blockNo)];
%     dlmwrite(fil_nm,temp,'delimiter','','-append')
%     temp=['A' num2str(blockNo) ', C' num2str(blockNo) ', ' num2str(setABCD(1)) ', ' num2str(setABCD(2))];
%     dlmwrite(fil_nm,temp,'precision',12,'delimiter','','-append')
%     
%     temp=['*Nfill, NSET=EG',num2str(blockNo)];
%     dlmwrite(fil_nm,temp,'delimiter','','-append')
%     temp=['E' num2str(blockNo) ', G' num2str(blockNo) ', ' num2str(setABCD(1)) ', ' num2str(setABCD(2))];
%     dlmwrite(fil_nm,temp,'precision',12,'delimiter','','-append')
%     
%     temp=['*Nfill, NSET=ACEG',num2str(blockNo)];
%     dlmwrite(fil_nm,temp,'delimiter','','-append')
%     temp=['AC' num2str(blockNo) ', EG' num2str(blockNo) ', ' num2str(setABCDEFGH(1)) ', ' num2str(setABCDEFGH(2))];
%     dlmwrite(fil_nm,temp,'precision',12,'delimiter','','-append')
% end

% Define the top node base surface for contact

%     temp=['***********************************************'];
%     dlmwrite(fil_nm,temp,'delimiter','','-append')
%     temp=['** Top Face Node Set'];
%     dlmwrite(fil_nm,temp,'delimiter','','-append')
%     
%     temp=['*Nfill, NSET=BD',num2str(blockNo)];
%     dlmwrite(fil_nm,temp,'delimiter','','-append')
%     temp=['B' num2str(blockNo) ', D' num2str(blockNo) ', ' num2str(setABCD(1)) ', ' num2str(setABCD(2))];
%     dlmwrite(fil_nm,temp,'precision',12,'delimiter','','-append')
%     
%     temp=['*Nfill, NSET=FH',num2str(blockNo)];
%     dlmwrite(fil_nm,temp,'delimiter','','-append')
%     temp=['F' num2str(blockNo) ', H' num2str(blockNo) ', ' num2str(setABCD(1)) ', ' num2str(setABCD(2))];
%     dlmwrite(fil_nm,temp,'precision',12,'delimiter','','-append')
%     
%     temp=['*Nfill, NSET=BDFH',num2str(blockNo)];
%     dlmwrite(fil_nm,temp,'delimiter','','-append')
%     temp=['BD' num2str(blockNo) ', FH' num2str(blockNo) ', ' num2str(setABCDEFGH(1)) ', ' num2str(setABCDEFGH(2))];
%     dlmwrite(fil_nm,temp,'precision',12,'delimiter','','-append')
    

%% 3] define elements
temp=['***********************************************'];
dlmwrite(fil_nm,temp,'delimiter','','-append')
temp=['*ELEMENT,TYPE=C3D8'];
dlmwrite(fil_nm,temp,'delimiter','','-append')
dlmwrite(fil_nm,Element1,'precision',12,'delimiter',',','-append')

temp=['*ELGEN, elset=allEL'];
dlmwrite(fil_nm,temp,'delimiter','','-append')
dlmwrite(fil_nm,ELGEN,'precision',12,'delimiter',',','-append')

% % Define top element base surface
% topEl = Element1(1) + ELGEN(7)*(ELGEN(5)-1);
% topElNd = Element1(2:end) + ((ELGEN(5)-1)*ELGEN(6));
% TpElNd = [topEl topElNd];
% 
% surfELs = [topEl ELGEN(2:4) ELGEN(8:10)];
% 
% temp=['***********************************************'];
% dlmwrite(fil_nm,temp,'delimiter','','-append')
% temp=['*ELEMENT,TYPE=C3D8'];
% dlmwrite(fil_nm,temp,'delimiter','','-append')
% dlmwrite(fil_nm,TpElNd,'precision',12,'delimiter',',','-append')
% 
% temp=['*ELGEN, elset=TopSurfElm'];
% dlmwrite(fil_nm,temp,'delimiter','','-append')
% dlmwrite(fil_nm,surfELs,'precision',12,'delimiter',',','-append')



