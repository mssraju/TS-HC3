function PSize = SubregionSizes(p1, cluster, N)
% Allocate N solutions to the subregions in proportion to the hypervolume
% of each cluster. Clusters with fewer than 2 members (HV = 0) receive half
% of the smallest positive HV, or 0.5 if every cluster has HV = 0.

    K  = length(cluster);
    hv = zeros(1, K);
    for i = 1:K
        members = cluster(i).sol;
        if numel(members) >= 2
            hv(i) = HV(p1(members), p1.objs);   % exact for M < 4
        end
    end

    zeroIdx = find(hv == 0);
    posIdx  = setdiff(1:K, zeroIdx);
    minPos  = min(hv(posIdx));
    if length(zeroIdx) == K
        hv(zeroIdx) = 0.5;
    else
        hv(zeroIdx) = minPos / 2;
    end

    PSize = round(hv / sum(hv) * N);
end