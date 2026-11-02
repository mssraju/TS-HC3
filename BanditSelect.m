function action = BanditSelect(Q, numStrategies, gStart, maxGen, gen)
% Decaying epsilon-greedy arm selection.
% Random choice during the first half of the remaining generations, and with
% probability epsilon afterwards; otherwise the arm with the largest reward.

    epsilon0 = 0.1;      % initial exploration rate
    decay    = 0.001;    % decay rate
    epsilon  = epsilon0 / (1 + decay * gen);

    if gStart + (maxGen - gStart)/2 > gen || rand < epsilon
        action = randi(numStrategies);
    else
        [~, action] = max(Q);
    end
end