function Q = reward(Q, Off, action, retainedRatio)
% Add the reward of the chosen arm: the mean of
%   (i)  the feasible ratio of the offspring, and
%   (ii) the normalised number of offspring retained by the populations.
% The two terms are weighted equally (0.5 / 0.5).

    CV            = sum(max(0, Off.cons), 2);
    feasibleRatio = sum(CV == 0) / length(Off);
    Q(action)     = Q(action) + mean([feasibleRatio, retainedRatio]);
end