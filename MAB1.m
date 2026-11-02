function action = MAB1(Q, numStrategies, gStart, maxGen, gen)
% Level 2: choose the variation operator (arm 1 = DE_Feasible, arm 2 = GA/DE mix).
    action = BanditSelect(Q, numStrategies, gStart, maxGen, gen);
end