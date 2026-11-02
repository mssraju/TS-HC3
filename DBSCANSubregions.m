function [cluster, Center] = DBSCANSubregions(pp1, numObj)
% Cluster the front with DBSCAN (cosine distance) on normalised objectives.
%   minPts : max(2, ceil(0.02*|X|)) for M>=3, and 1 for M=2
%   eps    : chosen from a 200-point grid by the best average silhouette;
%            falls back to 0.1 if no grid value gives a valid clustering
%   Clusters are sorted by the last objective of their centre (descending).

    M = size(pp1, 2);
    [X, ~] = Normalization4(pp1, min(pp1), max(pp1));

    if numObj >= 3
        minPts = max(2, ceil(size(X,1) * 0.02));
    else
        minPts = 1;
    end

    epsGrid   = linspace(0.001, 1, 200);
    bestScore = -Inf;
    bestEps   = NaN;
    for i = 1:length(epsGrid)
        labels = dbscan(X, epsGrid(i), minPts, 'Distance', 'cosine');
        valid  = labels(labels > 0);
        if numel(unique(valid)) > 1 && ~ismember(-1, labels)
            try
                silh  = silhouette(X, labels, 'cosine');
                score = mean(silh(~isnan(silh)));
                if score > bestScore
                    bestScore = score;
                    bestEps   = epsGrid(i);
                end
            catch
                continue;
            end
        end
    end
    if isnan(bestEps)
        bestEps = 0.1;
    end

    idx = dbscan(X, bestEps, minPts, 'Distance', 'cosine');
    uc  = unique(idx(idx > 0));
    if isempty(uc)                      % no cluster found: treat everything as one
        uc  = 1;
        idx = ones(size(pp1,1), 1);
    end

    % Members and centres (centres in the original objective space)
    K       = length(uc);
    cluster = struct('sol', cell(1,K));
    Center  = zeros(K, M);
    for i = 1:K
        members        = find(idx == uc(i));
        cluster(i).sol = members;
        if length(members) > 1
            Center(i,:) = mean(pp1(members,:), 1);
        else
            Center(i,:) = pp1(members,:);
        end
    end

    [~, order] = sort(Center(:,end), 'descend');
    cluster    = cluster(order);
    Center     = Center(order,:);
end