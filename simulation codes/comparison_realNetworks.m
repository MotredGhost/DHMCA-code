clc
clear
%% Number of nodes and medoids, trials and the order.
totalTrials = 10000; 
performanceDistributed = zeros(1, totalTrials);
performancefastPAM = zeros(1, totalTrials);
performancePAM1 = zeros(1, totalTrials);
performanceOldPAM = zeros(1, totalTrials);
diameter = zeros(1, totalTrials);
numberOfTrials = 0;
for trial = 1 : totalTrials
    if numberOfTrials == 100
        break
    end
%% flag = 1: use power network, otherwise use distributed water network
load('power_networks1.mat');
flag = 1;
if flag == 1
ajacency_matrix = GD98_112nodes_168edges_3degrees_AdjM;
else
ajacency_matrix = readmatrix('distributed_water_networks');
end
n = length(ajacency_matrix(1,:));
the_graph = graph(ajacency_matrix);
Neighbors = zeros(n, n);
for i = 1 : n
len = length(neighbors(the_graph, i));
Neighbors(i, 1 : len) = neighbors(the_graph, i)';
end    
M = 7;  % number of medoids
K = 1;  % each cluster has a number of "K" medoids
round = 1000;
performance = zeros(1, round);
%performance1 = zeros(1, round);
d_true = zeros(n, K);   % record distance information
m_true = zeros(n, K);    % record the closest "K" medoids
p_true = zeros(n, K);    % record parent node for a medoid
medoidIndexArray = zeros(round + 1, M);  % record the medoids for each run

%% generate medoids
flag = 0;
while flag ~= 1
medoidIndexArray(1,:) = randi([1 n],1,M);
if  length(unique(medoidIndexArray(1,:))) == length(medoidIndexArray(1,:))
    flag = 1;   % all initial medoids are distinct
end
end
medoidStatus = zeros(1,n);
for i = 1 : M
    medoidStatus(1, medoidIndexArray(1,i)) = 1;
end
flag = graphConnectedOrNot(Neighbors);
if flag == 0
    continue
else
    numberOfTrials = numberOfTrials + 1;
    fprintf('the %d-th effective trial\n', numberOfTrials);
end
s = zeros(1, n^2);
tt = zeros(1, n^2);
count = 1;
for i = 1 : n
    for j = 1 : length(nonzeros(Neighbors(i,:)))
        tt(1, count) = i;
        s(1,count) = Neighbors(i,j);
        count = count + 1;
    end
    
end
tt = nonzeros(tt)';
s = nonzeros(s)';
G = graph(s,tt);
Dmatrix = distances(G);
[idx, C, sumd, D] = kmedoids( (1:n)', M, 'Distance', @(ZI, ZJ) Dmatrix(ZI,ZJ), 'Algorithm','small' );
performancefastPAM(1, trial) = sum(sumd);
[idx1,C1,TD,ite] = FastPAM1(Dmatrix,M);
performancePAM1(1, trial) = TD;
[idx2, C2, sumd2, D2] = kmedoids( (1:n)', M, 'Distance', @(ZI, ZJ) Dmatrix(ZI,ZJ));
performanceOldPAM(1, trial) = sum(sumd2);
diameter(1, trial) = max(max(Dmatrix));
for time = 1 : round
%% state variables for high-order shortest path algorithm
d = zeros(n, K, diameter(1, trial) + 2);    % distance estimate
m = zeros(n, K, diameter(1, trial) + 2);    % medoids estimate
p = zeros(n, K, diameter(1, trial) + 2);    % parent estimate

%% initial states
for i = 1 : n
    m(i,:,1) = i * ones(1, K);
    if medoidStatus(1, i) == 1
        d(i,1,1) = 0;
        d(i,2:K,1) = randi([1,n],1,K-2+1);
    else
        d(i,:,1) = randi([1,n],1,K);
    end
end
%% high-order shortest path algorithm
for iteration = 2 : diameter(1, trial) + 2
    for i = 1 : n  
                 if medoidStatus(1, i) == 1
                    d(i,1,iteration) = 0;
                    m(i,1,iteration) = i;
                    p(i,1,iteration) = i;
                 else
                    len = length( nonzeros(Neighbors(i,:)) );
                    temp = intmax;
                    temp1 = 0;
                    parentInfo = zeros(1,2);
                    %% first find the minimum distance among neighbors
                            for l1 = 1 : len
                                for l2 = 1 : K
                                        if d(Neighbors(i,l1), l2, iteration - 1) < temp
                                            parentInfo(1,1) = Neighbors(i,l1);
                                            parentInfo(1,2) = l2;
                                            temp = d(Neighbors(i,l1), l2, iteration - 1); 
                                            temp1 = m(Neighbors(i,l1), l2, iteration - 1);
                                        elseif d(Neighbors(i,l1), l2, iteration - 1) == temp && m(Neighbors(i,l1), l2, iteration - 1) > temp1
                                            parentInfo(1,1) = Neighbors(i,l1);
                                            parentInfo(1,2) = l2;
                                            temp = d(Neighbors(i,l1), l2, iteration - 1); 
                                            temp1 = m(Neighbors(i,l1), l2, iteration - 1);
                                        elseif d(Neighbors(i,l1), l2, iteration - 1) == temp && m(Neighbors(i,l1), l2, iteration - 1) == temp1 ...
                                                && ( Neighbors(i,l1) >  parentInfo(1,1) || (Neighbors(i,l1) == parentInfo(1,1) && l2 > parentInfo(1,2) ) )  
                                            parentInfo(1,1) = Neighbors(i,l1);
                                            parentInfo(1,2) = l2;
                                            temp = d(Neighbors(i,l1), l2, iteration - 1); 
                                            temp1 = m(Neighbors(i,l1), l2, iteration - 1);                                        
                                        end
                                end
                            end

                    %% update the distance estimate 
                        
                            d(i,1,iteration) = d(parentInfo(1,1), parentInfo(1,2), iteration - 1) + 1;
                            m(i,1,iteration) = m(parentInfo(1,1), parentInfo(1,2), iteration - 1);
                            p(i,1,iteration) = parentInfo(1,1);
                 end
        for k = 2 : K
                %% find candidate neighbors
                      %len = length( nonzeros(Neighbors(i,:)) );
                      temp = intmax;
                      temp1 = 0;
                      parentInfo1 = zeros(1,2);
                      Mik = [m(i, 1 : k - 1, iteration), i];
                            for l1 = 1 : len
                                for l2 = 1 : K
                                    if ismember(m(Neighbors(i,l1), l2, iteration - 1), Mik) == 0 
                                        if d(Neighbors(i,l1), l2, iteration - 1) < temp
                                            parentInfo1(1,1) = Neighbors(i,l1);
                                            parentInfo1(1,2) = l2;
                                            temp = d(Neighbors(i,l1), l2, iteration - 1); 
                                            temp1 = m(Neighbors(i,l1), l2, iteration - 1);
                                        elseif d(Neighbors(i,l1), l2, iteration - 1) == temp && m(Neighbors(i,l1), l2, iteration - 1) > temp1
                                            parentInfo1(1,1) = Neighbors(i,l1);
                                            parentInfo1(1,2) = l2;
                                            temp = d(Neighbors(i,l1), l2, iteration - 1); 
                                            temp1 = m(Neighbors(i,l1), l2, iteration - 1);
                                        elseif d(Neighbors(i,l1), l2, iteration - 1) == temp && m(Neighbors(i,l1), l2, iteration - 1) == temp1 ...
                                                && ( Neighbors(i,l1) >  parentInfo1(1,1) || (Neighbors(i,l1) == parentInfo1(1,1) && l2 > parentInfo1(1,2) ) )  
                                            parentInfo1(1,1) = Neighbors(i,l1);
                                            parentInfo1(1,2) = l2;
                                            temp = d(Neighbors(i,l1), l2, iteration - 1); 
                                            temp1 = m(Neighbors(i,l1), l2, iteration - 1);                                        
                                        end
                                    end
                                end
                            end
                 %% update the distance estimate 
                 if parentInfo1(1,1) == 0    % if no such neighbor
                            d(i,k,iteration) = n;
                            m(i,k,iteration) = i;
                            p(i,k,iteration) = i;
                 else
                            d(i,k,iteration) = d(parentInfo1(1,1), parentInfo1(1,2), iteration - 1) + 1;
                            m(i,k,iteration) = m(parentInfo1(1,1), parentInfo1(1,2), iteration - 1);
                            p(i,k,iteration) = parentInfo1(1,1);
                 end
        end
    end
      
end
%% store true distance, closest medoids and parents
d_true(:,:) = d(:,:,  iteration );
m_true(:,:) = m(:,:,iteration );
p_true(:,:) = p(:,:,iteration );
%% state variables for the aggregation algorithm
O = zeros(n, K, diameter(1, trial) + 2); % second position refers to the corresponding medoid, e.g., O(3, 1, 5): the aggregation value of node 3 with respect to its closest medoid at the 5-th iteration
numDescendents = zeros(n, K, diameter(1, trial) + 2);
for iteration = 2 : diameter(1, trial) + 2
    for i = 1 : n
        len = length(nonzeros(Neighbors(i,:)));
        count = 0;
        Descendents = zeros(n, 3);             % descendent of node i; 1st element: descendent ID, 2nd element: order of the medoid that the neigboring node takes i as the parent; 3rd element: the corresponding medoid ID
                for l1 = 1:len                 % Descendents(3, 2, 260) means the node with index 3 takes i as the parent w.r.t. the medoid/node with index 260, and such medoid is its 2nd closest medoid
                    position = find(p_true(Neighbors(i,l1), :) == i);  % if this neighbor takes i as a parent; a neighboring may take i as the parent w.r.t. different medoids
                    if length(position) ~= 0
                        for l0 = 1 : length(position)
                        count = count + 1;
                        Descendents(count,:) = [Neighbors(i,l1) position(1, l0) m_true( Neighbors(i,l1), position(1, l0) ) ];
                        end                
                    end        
                end
                %% calculate number of descendents, sum of distances
                for l2 = 1 : K
                    for l3 = 1 : count
                        if Descendents(l3, 3) == m_true(i, l2)
                            O(i, l2, iteration) = O(i, l2, iteration) + O( Descendents(l3,1), Descendents(l3,2), iteration - 1 );
                            numDescendents(i, l2, iteration) = numDescendents(i, l2, iteration) + numDescendents( Descendents(l3,1), Descendents(l3,2), iteration - 1 );
                        end
                    end
                    O(i, l2, iteration) = O(i, l2, iteration) + d_true(i, l2);           % add i's own value 
                    numDescendents(i, l2, iteration) = numDescendents(i, l2, iteration) + 1;   % count i itself 
                end
    end  
end
%% performance function
performance(1,time) = performance(1,time) + sum(O(medoidIndexArray(time, :), 1, iteration));
if time >= 3 && performance(1,time) == performance(1,time - 1)
    break
end
%% find children of each medoid w.r.t. the medoid itself
medoidHasBeenSelected = zeros(1, M);
for i = 1 : M
    count = 0;
    flag = 0;    % 0: there is a new medoid added
    childrenInfo = zeros(n,2);
    theMedoid = medoidIndexArray(time, i);
    len = length(nonzeros(Neighbors(theMedoid,:)));
    for l1 = 1 : len
        position2 = find( m_true(Neighbors(theMedoid, l1), :) == theMedoid );
        len1 = length(position2);
        if len1 > 1
            fprintf('error!!!!!!!!!!')
        end
        if  len1 ~= 0
            count = count + 1;
            % the "numDescendents" here includes the node itself, so we need to decrease "numDescendents" by one  
            childrenInfo(count, :) = [Neighbors(theMedoid, l1)  numDescendents(Neighbors(theMedoid, l1), position2, iteration) - 1];  % child ID and the number of its descendents
        end
    end
    childrenInfo( all(~childrenInfo,2), : ) = [];
    %% update medoid
    count1 = 1;
    len1 = length(childrenInfo(:,1));      % number of children of a medoid
    if len1 ~= 0
    for l2 = 1 : len1
        if childrenInfo(l2,2) >= sum(childrenInfo(:,2)) - childrenInfo(l2,2) + len1  && ismember(childrenInfo(l2,1), medoidHasBeenSelected) ~= 1 ...
               && ismember(childrenInfo(l2,1), medoidIndexArray(time,:)) ~= 1   % 1. not has been selected as a medoid; 2. is not a current medoid
            medoidStatus(1, childrenInfo(l2,1)) = 1;
            medoidStatus(1, theMedoid) = 0;
            flag = 1;
            temp = childrenInfo(l2,1);
            break
        end
    end  
    end
            if flag == 1
            medoidHasBeenSelected(1, count1) = temp;
            count1 = count1 + 1;
            end
end
count = 1;
for i = 1 : n
    if  medoidStatus(1, i) == 1
        medoidIndexArray(time + 1, count) = i;
        count = count + 1;
    end
end

end
performanceDistributed(1, trial) = performance(1,time);
end
fprintf('the avreage result of distributed algorithm is %d\n', mean(nonzeros(performanceDistributed)));
fprintf('the average result of PAM with k-means++ intialization is %d\n', mean(nonzeros(performanceOldPAM)));
fprintf('the average result of fastPAM  is %d\n', mean(nonzeros(performancefastPAM)));
fprintf('the average result of newPAM is %d\n', mean(nonzeros(performancePAM1)));
fprintf('the ratio is %d\n', mean(nonzeros(performancePAM1))/mean(nonzeros(performanceDistributed)) );

%fprintf('the average diameter is %d\n', mean(nonzeros(diameter)));

