function [dist,path] = dijkstra(nodes,segments,start_id,finish_id)
%DIJKSTRA Calculates the shortest distance and path between points on a map
%   using Dijkstra's Shortest Path Algorithm
% 
% [DIST, PATH] = DIJKSTRA(NODES, SEGMENTS, SID, FID)
%   Calculates the shortest distance and path between start and finish nodes SID and FID
% 
% [DIST, PATH] = DIJKSTRA(NODES, SEGMENTS, SID)
%   Calculates the shortest distances and paths from the starting node SID to all
%     other nodes in the map
% 
% Note:
%     DIJKSTRA is set up so that an example is created if no inputs are provided,
%       but ignores the example and just processes the inputs if they are given.
% 
% Inputs:
%     NODES should be an Nx3 or Nx4 matrix with the format [ID X Y] or [ID X Y Z]
%       where ID is an integer, and X, Y, Z are cartesian position coordinates)
%     SEGMENTS should be an Mx3 matrix with the format [ID N1 N2]
%       where ID is an integer, and N1, N2 correspond to node IDs from NODES list
%       such that there is an [undirected] edge/segment between node N1 and node N2
%     SID should be an integer in the node ID list corresponding with the starting node
%     FID (optional) should be an integer in the node ID list corresponding with the finish
% 
% Outputs:
%     DIST is the shortest Euclidean distance
%       If FID was specified, DIST will be a 1x1 double representing the shortest
%         Euclidean distance between SID and FID along the map segments. DIST will have
%         a value of INF if there are no segments connecting SID and FID.
%       If FID was not specified, DIST will be a 1xN vector representing the shortest
%         Euclidean distance between SID and all other nodes on the map. DIST will have
%         a value of INF for any nodes that cannot be reached along segments of the map.
%     PATH is a list of nodes containing the shortest route
%       If FID was specified, PATH will be a 1xP vector of node IDs from SID to FID.
%         NAN will be returned if there are no segments connecting SID to FID.
%       If FID was not specified, PATH will be a 1xN cell of vectors representing the
%         shortest route from SID to all other nodes on the map. PATH will have a value
%         of NAN for any nodes that cannot be reached along the segments of the map.
% 
% Example:
%     dijkstra; % calculates shortest path and distance between two nodes
%               % on a map of randomly generated nodes and segments
% 
% Example:
%     nodes = [(1:10); 100*rand(2,10)]';
%     segments = [(1:17); floor(1:0.5:9); ceil(2:0.5:10)]';
%     figure; plot(nodes(:,2), nodes(:,3),'k.');
%     hold on;
%     for s = 1:17
%         if (s <= 10) text(nodes(s,2),nodes(s,3),[' ' num2str(s)]); end
%         plot(nodes(segments(s,2:3)',2),nodes(segments(s,2:3)',3),'k');
%     end
%     [d, p] = dijkstra(nodes, segments, 1, 10)
%     for n = 2:length(p)
%         plot(nodes(p(n-1:n),2),nodes(p(n-1:n),3),'r-.','linewidth',2);
%     end
%     hold off;
% 
% Author: Joseph Kirk
% Email: jdkirk630 at gmail dot com
% Release: 1.3
% Release Date: 5/18/07

if (nargin < 3) % SETUP
    % (GENERATE RANDOM EXAMPLE OF NODES AND SEGMENTS IF NOT GIVEN AS INPUTS)
    % Create a random set of nodes/vertices,and connect some of them with
    % edges/segments. Then graph the resulting map.
    num_nodes = 40; L = 100; max_seg_length = 30; ids = (1:num_nodes)';
    nodes = [ids L*rand(num_nodes,2)]; % create random nodes
    h = figure; plot(nodes(:,2),nodes(:,3),'k.') % plot the nodes
    text(nodes(num_nodes,2),nodes(num_nodes,3),...
        [' ' num2str(ids(num_nodes))],'Color','b','FontWeight','b')
    hold on
    num_segs = 0; segments = zeros(num_nodes*(num_nodes-1)/2,3);
    for i = 1:num_nodes-1 % create edges between some of the nodes
        text(nodes(i,2),nodes(i,3),[' ' num2str(ids(i))],'Color','b','FontWeight','b')
        for j = i+1:num_nodes
            d = sqrt(sum((nodes(i,2:3) - nodes(j,2:3)).^2));
            if and(d < max_seg_length,rand < 0.6)
                plot([nodes(i,2) nodes(j,2)],[nodes(i,3) nodes(j,3)],'k.-')
                % add this link to the segments list
                num_segs = num_segs + 1;
                segments(num_segs,:) = [num_segs nodes(i,1) nodes(j,1)];
            end
        end
    end
    segments(num_segs+1:num_nodes*(num_nodes-1)/2,:) = [];
    axis([0 L 0 L])
    % Calculate Shortest Path Using Dijkstra's Algorithm
    % Get random starting/ending nodes,compute the shortest distance and path.
    start_id = ceil(num_nodes*rand); disp(['start id = ' num2str(start_id)]);
    finish_id = ceil(num_nodes*rand); disp(['finish id = ' num2str(finish_id)]);
    [distance,path] = dijkstra(nodes,segments,start_id,finish_id);
    disp(['distance = ' num2str(distance)]); disp(['path = [' num2str(path) ']']);
    % If a Shortest Path exists,Plot it on the Map.
    figure(h)
    for k = 2:length(path)
        m = find(nodes(:,1) == path(k-1));
        n = find(nodes(:,1) == path(k));
        plot([nodes(m,2) nodes(n,2)],[nodes(m,3) nodes(n,3)],'ro-','LineWidth',2);
    end
    title(['Shortest Distance from ' num2str(start_id) ' to ' ...
        num2str(finish_id) ' = ' num2str(distance)])
    hold off
    
else %--------------------------------------------------------------------------
    % MAIN FUNCTION - DIJKSTRA'S ALGORITHM
    
    % initializations
    node_ids = nodes(:,1);
    [num_map_pts,cols] = size(nodes);
    table = sparse(num_map_pts,2);
    shortest_distance = Inf(num_map_pts,1);
    settled = zeros(num_map_pts,1);
    path = num2cell(NaN(num_map_pts,1));
    col = 2;
    pidx = find(start_id == node_ids);
    shortest_distance(pidx) = 0;
    table(pidx,col) = 0;
    settled(pidx) = 1;
    path(pidx) = {start_id};
    if (nargin < 4) % compute shortest path for all nodes
        while_cmd = 'sum(~settled) > 0';
    else % terminate algorithm early
        while_cmd = 'settled(zz) == 0';
        zz = find(finish_id == node_ids);
    end
    while eval(while_cmd)
        % update the table
        table(:,col-1) = table(:,col);
        table(pidx,col) = 0;
        % find neighboring nodes in the segments list
        neighbor_ids = [segments(node_ids(pidx) == segments(:,2),3);
            segments(node_ids(pidx) == segments(:,3),2)];
        % calculate the distances to the neighboring nodes and keep track of the paths
        for k = 1:length(neighbor_ids)
            cidx = find(neighbor_ids(k) == node_ids);
            if ~settled(cidx)
                d = sqrt(sum((nodes(pidx,2:cols) - nodes(cidx,2:cols)).^2));
                if (table(cidx,col-1) == 0) || ...
                        (table(cidx,col-1) > (table(pidx,col-1) + d))
                    table(cidx,col) = table(pidx,col-1) + d;
                    tmp_path = path(pidx);
                    path(cidx) = {[tmp_path{1} neighbor_ids(k)]};
                else
                    table(cidx,col) = table(cidx,col-1);
                end
            end
        end
        % find the minimum non-zero value in the table and save it
        nidx = find(table(:,col));
        ndx = find(table(nidx,col) == min(table(nidx,col)));
        if isempty(ndx)
            break
        else
            pidx = nidx(ndx(1));
            shortest_distance(pidx) = table(pidx,col);
            settled(pidx) = 1;
        end
    end
    if (nargin < 4) % return the distance and path arrays for all of the nodes
        dist = shortest_distance';
        path = path';
    else % return the distance and path for the ending node
        dist = shortest_distance(zz);
        path = path(zz);
        path = path{1};
    end
end