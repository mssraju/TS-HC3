function [wa,wb,znadTPop,znadArch,TPop,radius,PSize,numClusters,Center,r_circle,Gene,Population1] = ...
    Strat(Global,TPop,arch,znadTPop,znadArch,wa,wb,gen,radius,PSize,numClusters,Center,r_circle,Population1)
% Stage control and subregion construction.
%   wa = 1 : nadir point of the unconstrained population has stabilised
%   wb = 1 : nadir point of the feasibility population has stabilised ->
%            switch to stage 2 and build the first set of subregions
% Both tests compare the current nadir point with the one 'lookback'
% generations earlier.

    lookback = 50;      % generations between the two nadir points compared
    tol      = 1e-4;    % nadir-point change below which a population counts as stable

    znadTPop(gen,:) = max(TPop.objs, [], 1);
    Gene = 1;           % reset on every call; set to gen only when wa is activated

    %% Phase A: unconstrained population has converged (or half of the budget is used)
    if gen > lookback && wa == 0
        stable = pdist2(znadTPop(gen,:), znadTPop(gen-lookback,:), 'euclidean') <= tol;
        if stable || Global.maxFE/Global.FE <= 2
            wa   = 1;
            Gene = gen;
        end
    end

    %% Phase B: feasibility population has converged (or 55% of the budget is used)
    if wa == 1 && wb == 0
        if isempty(arch)
            znadArch(gen,:) = zeros(1, Global.M);
        else
            znadArch(gen,:) = max(arch.objs, [], 1);
        end

        stable = pdist2(znadArch(gen,:), znadArch(gen-lookback,:), 'euclidean') <= tol && length(arch) > 10;
        if stable || Global.FE/Global.maxFE >= 0.55
            wb          = 1;
            Population1 = TPop;     % unconstrained population becomes the epsilon-CHT population

            % Unique solutions of the feasibility population
            [~, id] = unique(arch.objs, 'rows');
            p1      = arch(id);
            pp1     = p1.objs;

            % Subregions: DBSCAN clusters, HV-based sizes, circular boundaries
            [cluster, Center] = DBSCANSubregions(pp1, Global.M);
            numClusters       = length(cluster);
            PSize             = SubregionSizes(p1, cluster, Global.N);
            [radius, r_circle] = cluster_radius(pp1, cluster, Center);

            % Initial subregion population
            inIdx = SubregionMembership(pp1, Center, r_circle);
            TPop  = SubregionSelection(p1, inIdx, PSize, @SPEA);
        end
    end
end