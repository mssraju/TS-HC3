function [Population, action] = MAB(Population1, Population2, Q, numStrategies, gStart, maxGen, gen)
% Level 1: choose the parent population (arm 1 = Population1, arm 2 = Population2).
    action = BanditSelect(Q, numStrategies, gStart, maxGen, gen);
    if action == 1
        Population = Population1;
    else
        Population = Population2;
    end
end