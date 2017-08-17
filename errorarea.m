function errorarea(x,y,l,u,col,face)
% Similar to errobar except produces shaded area around line for error
% function errorarea(x,y,l,u,col,face)
%
% Uses Matlab's plot function
%
% Steve Fleming 2012

if nargin < 5
    col = 'k';
    face = [0.6 0.6 0.6];
end

hold on
% Add error area as patch
for p = 1:length(x)-1
   
    pX = [x(p) x(p+1) x(p+1) x(p)];
    pY = [l(p) l(p+1) u(p+1) u(p)];
    h = fill(pX,pY,face);
    set(h,'edgecolor',face);
end

plot(x,y,'Color',col,'LineWidth',1.5);

plot(x,u,'Color',col,'LineStyle','--');
plot(x,l,'Color',col,'LineStyle','--');

hold off