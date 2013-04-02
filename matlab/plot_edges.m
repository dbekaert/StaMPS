function l=plot_edges(edges,x,y,z)
%PLOT_EDGES plot edges
%   function l=plot_edges(edges,x,y)
%   EDGES is Nx2 index to points in X and Y
%   X is Mx1 vector of x coords
%   Y is Mx1 vector of y coords
%   Z is Nx1 vector of attribute values or a colour (default 'b')
%
%   Andy Hooper June 2006

if nargin<4
    z='b';
end

if size(edges,2)>2
    edges=edges(:,end-1:end);
end

x1=x(edges(:,1));
x2=x(edges(:,2));
y1=y(edges(:,1));
y2=y(edges(:,2));

c=colormap;
plot(x,y,'r.')
l=line([x1';x2'],[y1';y2']);
if ~ischar(z)
    z=round(((z-min(z))/(max(z)-eps-min(z))*(size(c,1)-1))+0.5);
    for i=1:length(l)
        set(l(i),'color',c(z(i),:))
    end
else
    set(l,'color',z(1))
end
