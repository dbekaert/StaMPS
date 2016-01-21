function []=uw_interp();
%UW_INTERP Interpolate grid using nearest neighbour
%
%   Andy Hooper May 2007
%
%   ============================================================================
%   01/2012 AH: Speed up read/write for triangle 
%   01/2013 AH: Replace dsearch by dsearchn only for versions 2012 onwards
%   08/2014 DB: Suppress command line output
%   09/2015 AH: use matlab triangulation if triangle program not installed
%   ============================================================================

fprintf('Interpolating grid...\n')

uw=load('uw_grid','n_ps','n_ifg','nzix');

arch=computer('arch');
if strcmpi(arch(1:3),'win')
    use_triangle='n';
else
    tripath=system('which triangle >& /dev/null');
    if tripath==0
        use_triangle='y';
    else
        use_triangle='n';
    end  
end
    
[y,x]=find(uw.nzix);
xy=[[1:uw.n_ps]',x,y];

if use_triangle=='y'
    nodename=['unwrap.1.node'];
    fid=fopen(nodename,'w');
    fprintf(fid,'%d 2 0 0\n',uw.n_ps);

    fprintf(fid,'%d %d %d\n',xy');
    fclose(fid);

    [a,b] = system('triangle -e unwrap.1.node > triangle.log');

    fid=fopen('unwrap.2.edge','r');
    header=str2num(fgetl(fid));
    N=header(1);
    edgs=fscanf(fid,'%d %d %d %d\n',[4,N])';
    fclose(fid);
    n_edge=size(edgs,1);
    if n_edge~=N
        error('missing lines in unwrap.2.edge')
    end

    fid=fopen('unwrap.2.ele','r');
    header=str2num(fgetl(fid));
    N=header(1);
    ele=fscanf(fid,'%d %d %d %d\n',[4,N])';
    fclose(fid);
    n_ele=size(ele,1);
    if n_ele~=N
        error('missing lines in unwrap.2.ele')
    end
else
    xy=double(xy);
    ele=delaunay(xy(:,2),xy(:,3));
    tr=triangulation(ele,xy(:,2),xy(:,3));
    edgs=edges(tr);
    n_edge=size(edgs,1);
    edgs=[[1:n_edge]',edgs];
    n_ele=size(ele,1);
    ele=[[1:n_ele]',ele];
end    

z=[1:uw.n_ps];
[nrow,ncol]=size(uw.nzix);

[X,Y]=meshgrid(1:ncol,1:nrow);
matlab_version=version('-release');
if str2num(matlab_version(1:4))<2012
    Z=dsearch(x,y,ele(:,2:4),X,Y); % dsearch removed in MatlabR2012a
else
    Z=dsearchn([x,y],ele(:,2:4),[X(:),Y(:)]); %index from grid to pixel node
    Z = reshape(Z,nrow,ncol);
end
Zvec=Z(:);
grid_edges=[Zvec(1:end-nrow),Zvec(nrow+1:end)]; % col edges
Zvec=reshape(Z',nrow*ncol,1);
grid_edges=[grid_edges;[Zvec(1:end-ncol),Zvec(ncol+1:end)]]; % add row edges
[sort_edges,I_sort]=sort(grid_edges,2); % sort each edge to have lowest pixel node first
edge_sign=I_sort(:,2)-I_sort(:,1);
[alledges,I,J]=unique(sort_edges,'rows'); % grid_edges=alledges(J)
sameix=(alledges(:,1)==alledges(:,2));
alledges(sameix,:)=0; % set edges connecting identical nodes to (0,0)
[edgs,I2,J2]=unique(alledges,'rows');
n_edge=size(edgs,1)-1;
edgs=[[1:n_edge]',edgs(2:end,:)]; % drop (0,0)
gridedgeix=(J2(J)-1).*edge_sign; % index to edges
colix=reshape(gridedgeix(1:nrow*(ncol-1)),nrow,ncol-1);
rowix=reshape(gridedgeix(nrow*(ncol-1)+1:end),ncol,nrow-1)';

fprintf('   Number of unique edges in grid: %d\n',n_edge);


save('uw_interp','edgs','n_edge','rowix','colix','Z');
