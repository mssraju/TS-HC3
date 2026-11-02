function [Population,Z,count] = EPCHT(Population,Offspring,B,W,Problem,epsilon_k,Z,nr)
% Update of the epsilon-CHT population (MOEA/D-style neighbourhood replacement).
% Offspring i may replace up to nr solutions among its mating neighbours.
%   count : number of offspring that replaced at least one solution

    count = 0;
    Z = min(Z, min(Offspring.objs, [], 1));

    for i = 1:Problem.N
        % Neighbourhood (90%) or whole population (10%)
        if rand < 0.9
            P = B(i, randperm(size(B,2)));
        else
            P = randperm(Problem.N);
        end

        % Tchebycheff values and constraint violations
        gOld  = max(abs(Population(P).objs - repmat(Z, length(P), 1)) .* W(P,:), [], 2);
        gNew  = max(repmat(abs(Offspring(i).obj - Z), length(P), 1) .* W(P,:), [], 2);
        cvOld = overall_cv(Population(P).cons);
        cvNew = overall_cv(Offspring(i).con) * ones(length(P), 1);

        % Replacement rule
        bothWithinEps = (cvOld <= epsilon_k) & (cvNew <= epsilon_k);
        betterOrEqual = (gOld >= gNew) & (bothWithinEps | (cvOld == cvNew));
        replace       = find(betterOrEqual | (cvNew < cvOld), nr);

        if ~isempty(replace)
            count = count + 1;
            Population(P(replace)) = Offspring(i);
        end
    end
end