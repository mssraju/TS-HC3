function [SubPop,Center,radius,PSize,numClusters,r_circle,count] = EnvironmentalSelection_Cluster( ...
    Population,Offspring,arch,Problem,Center,radius,PSize,Gen,numClusters,r_circle)
% Subregion-aware environmental selection.
%   Population : cell array, one entry per subregion (previous generation)
%   Offspring  : offspring of the current generation
%   arch       : feasibility-driven population
%   Every 50 generations (and before 95% of the budget) the subregions are
%   rebuilt; otherwise the existing boundaries are reused.
%   count      : number of selected solutions that come from the offspring part

    rebuild = mod(Gen,50) == 0 && Problem.FE/Problem.maxFE <= 0.95;

    if rebuild
        % ---- Rebuild subregions from the non-dominated unique solutions ----
        if iscell(Population)
            oldSub = [Population{:}];
            merged = [arch, oldSub, Offspring];
        else
            oldSub = [];
            merged = arch;
        end

        [~, id] = unique(merged.objs, 'rows');
        p1      = merged(id);
        FrontNo = NDSort(p1.objs, p1.cons, 1);
        p1      = p1(FrontNo == 1);
        pp1     = p1.objs;

        [cluster, Center] = DBSCANSubregions(pp1, Problem.M);
        [radius, r_circle] = cluster_radius(pp1, cluster, Center);
        numClusters = length(cluster);
        PSize       = SubregionSizes(p1, cluster, Problem.N);

        pool     = [oldSub, p1, Offspring];
        numOld   = length([oldSub, p1]);
    else
        % ---- Reuse the existing boundaries ----
        oldSub = [Population{:}];
        pool   = [oldSub, arch, Offspring];
        numOld = length([oldSub, arch]);
    end

    % Assign the pool to subregions and select within each subregion
    inIdx    = SubregionMembership(pool.objs, Center, r_circle);
    selected = unique(vertcat(inIdx{:}));
    count    = sum(selected > numOld);
    SubPop   = SubregionSelection(pool, inIdx, PSize, @archive2REA);
end