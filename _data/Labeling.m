
L = bwlabeln(imerode(I==2,ones(4)));
max(L(:))
%%


mode(L(L~=0)) 

%%

[x,y,z] = ind2sub( size(L), find( L == 15));
plot3(x,y,z,'o')
axis equal
grid on
figure(gcf)