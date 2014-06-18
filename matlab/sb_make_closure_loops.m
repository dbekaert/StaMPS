function []=sb_make_closure_loops()
%SB_MAKE_CLOSURE_LOOPS
%   Produce a matrix containign 1,-1 or 0 elements definining which
%   interferograms are used in a closed loop test.
%   Ekbal Hussain @ Uni Leeds
%   25/10/2013
%
%   ======================================================================
%   06/2014 AH: Modify for StaMPS
%   ======================================================================
logit;
fprintf('Making closure loops...\n')

load psver
psname=['ps',num2str(psver)];
loopname=['phuw_loops',num2str(psver)];

ps=load(psname);

%Interferogram indexes
x = 1:ps.n_ifg; x = x';
intfg = [ps.ifgday_ix,x];          %where: col1=date 1, col2=dat2, col3=int_no.
clear x

%Remove interferograms dropped in previous stamps steps
drop_ifg = getparm('drop_ifg_index',1);
if ~isempty(drop_ifg)
    for n = 1:size(drop_ifg,2)
       intfg = intfg(find(intfg(:,3)~=drop_ifg(n)),:);
    end   
end

%% Create interferogram triangulalation matrix

%Acquisition date (node) indexes
iz = unique(intfg(:,1));

%Interferograms used
ints_used = [];
closed_loops = [];
for m=1:size(iz,1)
   
    %set node ID
    node_id = iz(m);

    iy = intfg(intfg(:,1)==node_id,:);      %Interferograms which connect this node 
    a = unique(iy(:,1:2));                  %The nodes used in the combinations

    if length(a) >= 3 

        %Create all possible combinations using these nodes and output the subset that
        %only begins with the selected start node.
        combinations = combnk(a,3);
        usecom = combinations((combinations(:,1)==node_id),:);

        int_list1 = [usecom(:,1) usecom(:,2)];  %Positive interfeograms in the loop
        int_list2 = [usecom(:,2) usecom(:,3)];  %Positive interfeograms in the loop
        int_list3 = [usecom(:,1) usecom(:,3)];  %Negative interfeograms in the loop

        %Find the interferogram index connecting the nodes
        %First positive interferogram
        temp_column = nan(size(int_list1,1),1);
        for k=1:size(int_list1,1)
           [~,~,ix] = intersect(int_list1(k,:),intfg(:,1:2),'rows');
           if isempty(ix)
               temp_column(k,1)=NaN;
           else
           temp_column(k,1)=intfg(ix,3);
           end
        end
        int_list1 = [int_list1 temp_column];
        inan1 = ~isnan(int_list1(:,3));
        clear temp_column temp

        %Second positive interferogram
        temp_column = nan(size(int_list2,1),1);
        for k=1:size(int_list2,1)
           [~,~,ix] = intersect(int_list2(k,:),intfg(:,1:2),'rows');
           if isempty(ix)
               temp_column(k,1)=NaN;
           else
           temp_column(k,1)=intfg(ix,3);
           end
        end
        int_list2 = [int_list2 temp_column];
        inan2 = ~isnan(int_list2(:,3));
        clear temp_column temp

        %The negative interferogram
        temp_column = nan(size(int_list3,1),1);
        for k=1:size(int_list3,1)
           [~,~,ix] = intersect(int_list3(k,:),intfg(:,1:2),'rows');
           temp_column(k,1)=intfg(ix,3);
        end
        int_list3 = [int_list3 temp_column];
        inan3 = ~isnan(int_list3(:,3));
        clear temp_column temp

        %Remove any rows that contain NaNs
        igood = [inan1, inan2, inan3];
        for k=1:size(igood,1)
            if igood(k,:) == [1 1 1];
                ixgood(k,1) = 1;
            else
                ixgood(k,1) = NaN;
            end
        end
        int_list1 = int_list1(~isnan(ixgood),:);
        int_list2 = int_list2(~isnan(ixgood),:);
        int_list3 = int_list3(~isnan(ixgood),:);
        clear inan1 inan2 inan3 ixgood

        %Create a vector of the interferograms used
        ints_used = unique([ints_used;int_list1(:,3);int_list2(:,3);int_list3(:,3);]);
        
        %Add these to a loop closure matrix
        closed_loops = [closed_loops; int_list1(:,3) int_list2(:,3) -int_list3(:,3);];
        clear int_list1 int_list2 int_list3
    end
end

clear combinations usecom iy a iz k igood

%% Higher order loops for those not fitted by a triangle

A = [closed_loops(:,1:2),-closed_loops(:,3)];
missing = intfg(~ismember(intfg(:,3),A),3);
clear A

if ~isempty(missing)
    
    %%% replace this code - use bperp from ps2
    load bperp.1.in
    bp = zeros(ps.n_image,1);
    bp(1:(ps.master_ix-1),1)=bperp_1(1:(ps.master_ix-1),1);   %Master bperp = 0
    bp((ps.master_ix+1):end,1)=bperp_1(ps.master_ix:end,1);
    %%%%

    Nodes = [(1:size(bp,1))', ps.day, bp];
    segments = [intfg(:,3), intfg(:,1:2)];

    %Remove the shortest link between the two nodes
    Segments = segments((segments(:,1)~=[missing]), :);

    higherO_loops = [];
    for m=1:size(missing,1)

        A = intfg(find(intfg(:,3)==missing(m,1)),1:2);

        %Calculate the shortest path that links the two nodes
        [~,path] = dijkstra(Nodes,Segments,A(:,1),A(:,2));

        %Find the interferogram index for this path
        for k=1:(size(path,2)-1)        
            if path(1,k+1) > path(1,k);
                temp_vector = [path(1,k),path(1,k+1)];
                [~,~,ix] = intersect(temp_vector,intfg(:,1:2),'rows');
                higherO_loops(m,k)=intfg(ix,3);
            else
                temp_vector = [path(1,k+1),path(1,k)];             
                [~,~,ix] = intersect(temp_vector,intfg(:,1:2),'rows');
                higherO_loops(m,k)=-intfg(ix,3);          
            end
        end

        %Add the final interfeogram to close the loop
        higherO_loops(m,size(path,2)) = -missing;
    end
    clear m A path k temp_vector ix Nodes segments Segments
end

%% The final closed network

if ~isempty(missing)
    intfg_network = nan((size(closed_loops,1)+size(higherO_loops,1)),size(higherO_loops,2));
   
    %Create the new matrix here. Add the triangular loops
    intfg_network(1:size(closed_loops,1),1:size(closed_loops,2)) = closed_loops;
    %Add the higher order loops
    intfg_network((size(closed_loops,1)+1):end,:) = higherO_loops;    
else
    intfg_network = closed_loops;
end

%Turn this into a full matrix including flags, where 1 = use and 0 = don't use
intfg_loops = zeros(size(intfg_network,1),size(intfg,1));

% This is not working 
for m=1:size(intfg_network,1);       
    for n=1:nnz(~isnan(intfg_network(m,:)));
        if intfg_network(m,n)<0;
            intfg_loops(m,(-1*intfg_network(m,n)))=-1;
        else
            intfg_loops(m,intfg_network(m,n))=1;
        end
    end  
end
logit(1);
save(loopname,'intfg_loops')





