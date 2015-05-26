function []=sb_make_closure_loops()
%SB_MAKE_CLOSURE_LOOPS
%   Produce a matrix containign 1,-1 or 0 elements definining which
%   interferograms are used in a closed loop test.
%   Ekbal Hussain @ Uni Leeds
%   25/10/2013
%
%   ======================================================================
%   06/2014 AH: Modify for StaMPS
%   08/2014 EH: Fixed bug for determining higher order loops
%   08/2014 EH: Ignores interferograms not connected in a closed loop, bug fixes
%   01/2015 EH: User input required for individual ints (not in a loop) connecting two looped networks
%   01/2015 EH: Bug fixes
%   04/2015 EH: Checks if you are running small baselines and not single master. Speed optimisations
%   04/2015 EH: Removed user input setion for ints connecting two loops
%   04/2015 EH: Allows for up to 12 interferograms in a loop
%   ======================================================================

logit;
fprintf('Making closure loops...\n')

load psver
psname=['ps',num2str(psver)];
loopname=['phuw_loops',num2str(psver)];

ps=load(psname);

% Check you are running small baselines and not single master
small_baseline_flag = getparm('small_baseline_flag');
if strcmpi(small_baseline_flag,'n')
    error('Error: This method does not work for single master. Use small baselines instead.')
end

%Interferogram indices
x = 1:ps.n_ifg; x = x';
intfg = [ps.ifgday_ix,x];          %where: col1=date 1, col2=dat2, col3=int_no.

%Remove interferograms dropped in previous stamps steps
drop_ifg = getparm('drop_ifg_index',1);
if ~isempty(drop_ifg)
    for n = 1:size(drop_ifg,2)
       intfg = intfg(intfg(:,3)~=drop_ifg(n),:);
    end   
end

% loading the single master bperp
small_baseline_flag = getparm('small_baseline_flag');
if strcmpi(small_baseline_flag,'y')
    temp = load('../psver');
    psname_ps = ['..' filesep 'ps' num2str(temp.psver) '.mat'];
    load(psname_ps,'bperp')
else
   load(psname,'bperp') 
end

clear x n

%% Create interferogram triangulalation matrix

%Acquisition date (node) indexes
iz = unique(intfg(:,1));

%Interferograms used
ints_used = [];
closed_loops = [];
for m=1:size(iz,1)
   
    %set node ID
    node_id = iz(m);

    iy = intfg(intfg(:,1)==node_id,:);      %Interferograms that connect to this node [nID nID intID]
    a = unique(iy(:,1:2));                  %The nodes used in the combinations

    if length(a) >= 3                       %Can be used in one or more loops

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
           if isempty(ix)
               temp_column(k,1)=NaN;
           else
               temp_column(k,1)=intfg(ix,3);
           end
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
        closed_loops = [closed_loops; int_list1(:,3) int_list2(:,3) -int_list3(:,3);];      %This array contains only int IDs
        clear int_list1 int_list2 int_list3      
    end
end

clear combinations usecom iy a iz k igood ix node_id small_baseline_flag

%% Higher order loops for those not fitted by a triangle

%Find the ints that were not fitted by a triangle
A = [closed_loops(:,1:2),-closed_loops(:,3)];
missing_loops = intfg(~ismember(intfg(:,3),A),:);

%Ignore unconnected single interferograms
missing = [];
if ~isempty(missing_loops)
    for i=1:size(missing_loops,1)
        occurences_node1  = sum(sum(intfg(:,1:2)==missing_loops(i,1)),2);
        occurences_node2  = sum(sum(intfg(:,1:2)==missing_loops(i,2)),2);
        if occurences_node1 ~= 1 && occurences_node2 ~= 1
            missing = vertcat(missing,missing_loops(i,3));
        end
    end
end

remove_from_missing = [];
if ~isempty(missing)
    
    % loading the single master perpendicular baselines
    temp = load('../psver');
    psname_ps = ['..' filesep 'ps' num2str(temp.psver) '.mat'];
    sm_bp = load(psname_ps,'bperp');
    bp = sm_bp.bperp;
    
    Nodes = [(1:size(bp,1))', ps.day, bp];
    segments = [intfg(:,3), intfg(:,1:2)];

    % Allow for loops containing up to 8 interferograms
    higherO_loops = nan(size(missing,1),12);
    for m=1:size(missing,1)

        A = intfg(intfg(:,3)==missing(m,1),1:2);
        
        %Remove the shortest link between the two nodes      
        Segments = segments;
        Segments(Segments(:,1) == missing(m,1),:) = [];
        
        %Calculate the shortest path that links the two nodes using
        %dijkstra's algorithm
        [~,path] = dijkstra(Nodes,Segments,A(:,1),A(:,2));
        
        %Remove individual ints (not in a looped network) connecting two looped networks
        if isnan(path)
            disp(['Interferogram ',num2str(missing(m,:)),' is not connected to a loop.']);
            
            %remove this interferogram from the missing list
            remove_from_missing = [remove_from_missing;m];
        else
            %close the loop
            path = [path, path(1,1)];

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
        end
    end
 
    %Remove rows and colums containing just nans
    higherO_loops = higherO_loops(~all(isnan(higherO_loops),2),~all(isnan(higherO_loops),1));
    %replace nans with zero (makes the next bit easier)
    higherO_loops(isnan(higherO_loops)) = 0;
    
    %Remove duplicated loops
    x = sort(abs(higherO_loops),2);
    y =  unique(x,'rows');
    [~,ix] = ismember(y,x,'rows');
    
    higherO_loops = higherO_loops(ix,:);
    %replace any zeros with nans
    higherO_loops(higherO_loops==0) = nan;
end

clear A i occurences_node1 occurences_node2 missing_loops psname_ps temp sm_bp bp m path k temp_vector ix Nodes segments Segments x y

%% The final closed network

%Remove the single interferograms connecting two loops
missing(remove_from_missing,:) = [];

if ~isempty(missing)
    intfg_network = nan((size(closed_loops,1)+size(higherO_loops,1)),size(higherO_loops,2));
   
    %Create the new matrix here. Add the triangular loops
    intfg_network(1:size(closed_loops,1),1:size(closed_loops,2)) = closed_loops;
    %Add the higher order loops
    intfg_network((size(closed_loops,1)+1):end,:) = higherO_loops;    
else
    intfg_network = closed_loops;
end

%Turn this into a full matrix including flags, where 1,-1 = use and 0 = don't use
%intfg_loops: rows = number of closed loops, columns = total no. of ints in network
intfg_loops = zeros(size(intfg_network,1),size(intfg,1));
 
for m=1:size(intfg_network,1);       
    for n=1:nnz(~isnan(intfg_network(m,:)));
        if intfg_network(m,n)<0;
            intfg_loops(m,(abs(intfg_network(m,n))))=-1;
        else
            intfg_loops(m,intfg_network(m,n))=1;
        end
    end  
end
logit(1);
save(loopname,'intfg_loops')





