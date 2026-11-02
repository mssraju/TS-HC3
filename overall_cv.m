function result = overall_cv(cons)
% Overall constraint violation: sum of the positive violations per solution.
    cons(cons <= 0) = 0;
    result = sum(abs(cons), 2);
end