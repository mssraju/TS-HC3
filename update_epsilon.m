function result = update_epsilon(relaxSpeed, epsilon_k, epsilon_0, rf, targetRatio, gen, Tc, cp)
% epsilon-CHT relaxation threshold.
%   rf < targetRatio : shrink the current threshold geometrically
%   otherwise        : follow a polynomial schedule from epsilon_0 to 0 at gen = Tc
% (the caller passes the current maximum constraint violation as epsilon_0)

    if rf < targetRatio
        result = (1 - relaxSpeed) * epsilon_k;
    else
        result = epsilon_0 * ((1 - gen/Tc) ^ cp);
    end
end