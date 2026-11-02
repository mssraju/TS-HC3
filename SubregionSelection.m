function SubPop = SubregionSelection(pool, inIdx, PSize, truncateFcn)
% For each subregion keep all of its solutions, or truncate to PSize(k) with
% truncateFcn (e.g. @SPEA or @archive2REA) when there are more.

    K      = length(inIdx);
    SubPop = cell(1, K);
    for k = 1:K
        members = pool(inIdx{k});
        if isempty(members)
            SubPop{k} = [];
        elseif length(members) > PSize(k)
            SubPop{k} = truncateFcn(members, PSize(k));
        else
            SubPop{k} = members;
        end
    end
end