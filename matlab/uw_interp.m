function []=uw_interp();
%UW_INTERP Interpolate grid using nearest neighbour
%
%   Andy Hooper May 2007
%


fprintf('Interpolating grid...\n')

uw=load('uw_grid');

[n_ps,n_ifg]=size(uw.ph);

disp(sprintf('   Number of interferograms: %d',n_ifg))
disp(sprintf('   Number of points per ifg: %d',n_ps))

nodename=['unwrap.1.node'];
fid=fopen(nodename,'w');
fprintf(fid,'%d 2 0 0\n',n_ps);


[y,x]=find(uw.nzix);

for i=1:n_ps
    fprintf(fid,'%d %d %d\n',i,x(i),y(i));
end

fclose(fid);

!triangle -e unwrap.1.node

fid=fopen('unwrap.2.edge','r');
header=str2num(fgetl(fid));
N=header(1);
edges=zeros(N,4);
for i=1:N
    edges(i,:)=str2num(fgetl(fid));
end
fclose(fid);

fid=fopen('unwrap.2.ele','r');
header=str2num(fgetl(fid));
N=header(1);
ele=zeros(N,4);
for i=1:N
    ele(i,:)=str2num(fgetl(fid));
end
fclose(fid);

z=[1:n_ps];
[nrow,ncol]=size(uw.nzix);


[X,Y]=meshgrid(1:ncol,1:nrow);
Z=dsearch(x,y,ele(:,2:4),X,Y); %index from grid to pixel node
Zvec=Z(:);
grid_edges=[Zvec(1:end-nrow),Zvec(nrow+1:end)]; % col edges
Zvec=reshape(Z',nrow*ncol,1);
grid_edges=[grid_edges;[Zvec(1:end-ncol),Zvec(ncol+1:end)]]; % add row edges
[sort_edges,I_sort]=sort(grid_edges,2); % sort each edge to have lowest pixel node first
edge_sign=I_sort(:,2)-I_sort(:,1);
[alledges,I,J]=unique(sort_edges,'rows'); % grid_edges=alledges(J)
sameix=(alledges(:,1)==alledges(:,2));
alledges(sameix,:)=0; % set edges connecting identical nodes to (0,0)
[edges,I2,J2]=unique(alledges,'rows');
n_edge=size(edges,1)-1;
edges=[[1:n_edge]',edges(2:end,:)]; % drop (0,0)
gridedgeix=(J2(J)-1).*edge_sign; % index to edges
colix=reshape(gridedgeix(1:nrow*(ncol-1)),nrow,ncol-1);
rowix=reshape(gridedgeix(nrow*(ncol-1)+1:end),ncol,nrow-1)';

fprintf('Number of unique edges in grid: %d\n\n',n_edge);


save('uw_interp','edges','n_edge','rowix','colix','Z');
