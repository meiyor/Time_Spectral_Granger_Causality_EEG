function p_adj = holm_pvals(p)
    p = p(:);
    [ps, idx] = sort(p);
    m = numel(ps);

    p_adj_sorted = (m - (1:m)' + 1) .* ps;

    % enforce monotonicity
    for i = 2:m
        p_adj_sorted(i) = max(p_adj_sorted(i), p_adj_sorted(i-1));
    end

    p_adj_sorted = min(p_adj_sorted, 1);

    p_adj = nan(size(p));
    p_adj(idx) = p_adj_sorted;
end