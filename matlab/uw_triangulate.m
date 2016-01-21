function []=UW_triangulate(ph,xy,day);
%UW_TRIANGULATE triangulate data

disp('Triangulating...')


[n_ps,n_ifg]=size(ph);

disp(sprintf('   Number of interferograms: %d',n_ifg))
disp(sprintf('   Number of points per ifg: %d',n_ps))




xy=double(xy);
ele=delaunay(xy(:,2),xy(:,3));
tr=triangulation(ele,xy(:,2),xy(:,3));
edgs=edges(tr);

n_edge=size(edgs,1);
edgs=[[1:n_edge]',edgs];

n_ele=size(ele,1);
ele=[[1:n_ele]',ele];
eledix=[1:n_ele]';
% create sparse matrix with nodes as indices and egde id as value
max_edge=max(max(edgs));    
edgeix=sparse([edgs(:,2);max_edge],[edgs(:,3);max_edge],[edgs(:,1);0]); % make sure edgeix has dimensions max_edge by max_edge

for i=1:n_ele;
    % for each triangle, lookup edge id for each combination of
    % nodes (1,2;2,3;3,1). Set to -ve if (2,1 or 3,2 or 1,3).
    eledix(i,2)=edgeix(ele(i,2),ele(i,3))-edgeix(ele(i,3),ele(i,2));
    eledix(i,3)=edgeix(ele(i,3),ele(i,4))-edgeix(ele(i,4),ele(i,3));
    eledix(i,4)=edgeix(ele(i,4),ele(i,2))-edgeix(ele(i,2),ele(i,4));
end

keep_ele_ix=(eledix(:,2)~=0&eledix(:,3)~=0&eledix(:,4)~=0);
ele=ele(keep_ele_ix,:);
ele(:,1)=[1:size(ele,1)]';
eledix=eledix(keep_ele_ix,:);
eledix(:,1)=[1:size(eledix,1)]';

%========================================================================
% create inverse index (points from edge to triangle). Max of 2 triangles
% per edge, if only 1, then both triangle indices set to same.
%========================================================================

eeixele=repmat(eledix(:,1),1,3);
eeixedge=eledix(:,2:4);
eeixele=eeixele(:).*sign(eeixedge(:));
eeixedge=abs(eeixedge(:));
eeixall=sortrows([eeixedge,eeixele]);
[edelix,I,J]=unique(eeixall(:,1));
edelix(:,3)=eeixall(I,2);
I=I-(diff([0;I])-1);
edelix(:,2)=eeixall(I,2);

if size(edelix,1)<size(edgs,1)
    error('there are orphaned edges - increase max_outer_len')
end

edge_nodes=edgs(:,2:3);
edge_nodes=edge_nodes(:);
edge_edge=repmat([1:n_edge]',2,1);
[sort_node,sort_ix]=sort(edge_nodes);
sort_edge=edge_edge(sort_ix);

save_node=n_ps;
i2_save=1;
for i2=1:n_edge*2
    if sort_node(i2)~=save_node
        ndelix{save_node,1}=sort_edge(i2_save:i2-1);
        save_node=sort_node(i2);
        i2_save=i2;
    end
end
ndelix{save_node,1}=sort_edge(i2_save:i2); % save last node

%ph(ph~=0)=ph(ph~=0)./abs(ph(ph~=0)); % just to make sure phase is normalized


save('uw_phase','ph','xy','edgs','ele','eledix','edelix','ndelix','n_ps','n_ifg','n_ele','n_edge','day');
