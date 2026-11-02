function inIdx = SubregionMembership(pop, Center, r_circle)
% Indices of the rows of pop (objective values) inside each subregion circle.
% A solution is assigned to every circle that contains it; solutions outside
% all circles are not assigned.

    K      = size(Center, 1);
    R      = max(r_circle(:), 0);
    D      = pdist2(pop, Center);            % N x K distances to the centres
    inMask = D <= (R.' + 1e-9);              % N x K logical
    inIdx  = arrayfun(@(k) find(inMask(:,k)), (1:K).', 'UniformOutput', false);
end