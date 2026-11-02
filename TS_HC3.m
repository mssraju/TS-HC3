classdef TS_HC3 < ALGORITHM
% <multi> <real/integer> <constrained>
% TS-HC3: Two-Stage Hybrid Tri-Population CMOEA with Subregion-Aware Selection
%
% Stage 1 (exploration):  an unconstrained population and a feasibility-driven
%                         population co-evolve (Strat decides when to stop).
% Stage 2 (exploitation): epsilon-CHT population, feasibility-driven population
%                         and a subregion population (DBSCAN-based) co-evolve.
%
%------------------------------- Reference --------------------------------
% Sri Srinivasa Raju Modampuri, et al., "A Two-Stage Evolutionary Algorithm
% with Subregion-Aware Selection for Constrained Multi-objective
% Optimization", Applied Soft Computing, 2026.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------
% NOTE: if your PlatEMO version expects the tag line "% <algorithm> <A>",
% keep that line from your original file instead of the one above.

    methods
        function main(Algorithm, Problem)
            %% Initialization
            % Population  : feasibility-driven population
            % TPop        : stage 1: unconstrained population
            %               stage 2: cell array, one entry per subregion
            % Population1 : stage 1: archive, stage 2: epsilon-CHT population
            [W, Problem, B, Population, Z, TPop, wa, wb, Population1, Z1, znadTPop, znadArch] = Initial(Problem);

            %% Parameters (set after Initial because Initial may change Problem.N)
            numStrategies = 2;      % arms in each MAB level
            relaxSpeed    = 0.05;   % epsilon-CHT: speed of reducing the relaxation
            cp            = 2;      % epsilon-CHT: reduction-rate exponent
            nr            = 2;      % epsilon-CHT: max replacements per offspring
            Tc            = 0.9 * ceil(Problem.maxFE / Problem.N);   % epsilon-CHT control generation
            maxGen        = Problem.maxFE / Problem.N;               % generation budget used by the MAB schedule

            %% State
            Gen         = 1;        % generation counter
            Gene        = 1;        % generation at which the exploitation stage starts (g_MAB)
            Center      = [];       % subregion centres
            radius      = [];       % subregion angular radii
            r_circle    = [];       % subregion circular radii
            PSize       = [];       % subregion sizes
            numClusters = [];
            count  = 0;             % offspring retained by the feasibility population
            count1 = 0;             % offspring retained by the epsilon-CHT population
            Q  = zeros(1, numStrategies);   % reward of each parent population
            Q1 = zeros(1, numStrategies);   % reward of each variation operator
            epsilon_k = max(overall_cv(Population.cons));

            %% Main loop
            while Algorithm.NotTerminated(Population)

                if wb == 0
                    % ---------------- Stage 1: exploration ----------------
                    Off       = AdapativeOffspring(TPop, Problem, B);
                    [TPop, ~] = Second_Stage_EnvironmentalSelection([TPop, Off], Problem.N, 2);
                    [Population, Z1, ~] = AssociationD([TPop, Population, Off], W, Z1, wb, Problem);

                    % Stage-transition test and subregion construction
                    [wa, wb, znadTPop, znadArch, TPop, radius, PSize, numClusters, Center, r_circle, Gene, Population1] = ...
                        Strat(Problem, TPop, Population, znadTPop, znadArch, wa, wb, Gen, ...
                              radius, PSize, numClusters, Center, r_circle, Population1);
                else
                    % ---------------- Stage 2: exploitation ----------------
                    % Offspring set 1: two-level MAB (parent population, then operator)
                    [Parent, action] = MAB(Population, Population1, Q, numStrategies, Gene, maxGen, Gen);
                    opIdx            = MAB1(Q1, numStrategies, Gene, maxGen, Gen);
                    Off1             = Adapative1(Parent, Problem, B, opIdx);

                    % Offspring set 2: DE_Feasible on the subregion population
                    Off2 = AdapativeOffspringSR1(MergeSubregions(TPop), Problem);
                    Off  = [Off1, Off2];

                    % Feasibility-driven population
                    [Population, Z1, count] = AssociationD([Population, Off], W, Z1, wb, Problem);

                    % epsilon-CHT population
                    cvPop1    = sum(max(0, Population1.cons), 2);
                    rf        = sum(cvPop1 <= 1e-6) / Problem.N;                       % feasible ratio
                    alphaT    = 2 ./ (1 + exp(-Problem.FE * 10 / Problem.maxFE)) - 1;  % target ratio
                    epsilon_k = update_epsilon(relaxSpeed, epsilon_k, max(cvPop1), rf, alphaT, Gen, Tc, cp);
                    [Population1, Z, count1] = EPCHT(Population1, Off, B, W, Problem, epsilon_k, Z, nr);

                    % Subregion-aware environmental selection
                    [TPop, Center, radius, PSize, numClusters, r_circle, count3] = ...
                        EnvironmentalSelection_Cluster(TPop, Off, Population, Problem, ...
                                                       Center, radius, PSize, Gen, numClusters, r_circle);

                    % MAB reward updates
                    retained = mean([count, count1, count3]) / Problem.N;
                    Q  = reward(Q,  Off, action, retained);
                    Q1 = reward(Q1, Off, opIdx,  retained);
                end

                % ---------------- Final population selection ----------------
                if Problem.FE >= Problem.maxFE
                    [Population, ~] = BestPop(Population, MergeSubregions(TPop), Population1);
                end

                Gen = Gen + 1;
            end
        end
    end
end

% Merge the per-subregion cell array into one population array
function P = MergeSubregions(TPop)
    if iscell(TPop)
        P = [TPop{:}];
    else
        P = TPop;
    end
end