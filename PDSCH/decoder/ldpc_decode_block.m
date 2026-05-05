function [decoded, success] = ldpc_decode_block(llr_in, max_iter, info_len)
    persistent cache
    
    llr_in = llr_in(:);
    
    if isempty(cache) || cache.info_len ~= info_len
        if info_len > 3824
            Kb = 22; nb = 68; mb = 46; bg = 1;
        else
            Kb = 10; nb = 52; mb = 42; bg = 2;
        end
        Z_set = sort([2 4 8 16 32 64 128 256 3 6 12 24 48 96 192 384 ...
                       5 10 20 40 80 160 320 7 14 28 56 112 224 ...
                       9 18 36 72 144 288 11 22 44 88 176 352 ...
                       13 26 52 104 208 15 30 60 120 240]);
        idx = find(Kb * Z_set >= info_len, 1);
        if isempty(idx), Z = Z_set(end); else, Z = Z_set(idx); end
        
        c.info_len = info_len;
        c.K_ldpc = Kb * Z;
        c.N = nb * Z;
        c.M = mb * Z;
        c.parity_len = c.N - c.K_ldpc;
        c.K_orig = length(llr_in) - c.parity_len;
        if c.K_orig < 1, c.K_orig = info_len; c.parity_len = length(llr_in) - c.K_orig; end
        
        c.H = build_H(nb, mb, Z, bg);
        [c.edge_r, c.edge_c] = find(c.H);
        c.ne = length(c.edge_r);
        
        % Pre-compute row groupings: for each check node, which edges belong to it
        % Sort edges by row for efficient vectorized processing
        [c.edge_r_sorted, sort_idx] = sort(c.edge_r);
        c.edge_c_sorted = c.edge_c(sort_idx);
        c.sort_idx = sort_idx;
        
        % Inverse sort: map from sorted position back to original
        c.inv_sort_idx = zeros(c.ne, 1);
        c.inv_sort_idx(sort_idx) = (1:c.ne)';
        
        % Row boundaries: row_start(mi) to row_end(mi) gives edges for check node mi
        c.row_start = zeros(c.M, 1);
        c.row_end = zeros(c.M, 1);
        cur = 1;
        for mi = 1:c.M
            c.row_start(mi) = cur;
            while cur <= c.ne && c.edge_r_sorted(cur) == mi
                cur = cur + 1;
            end
            c.row_end(mi) = cur - 1;
        end
        
        % Pre-compute column groupings similarly (sorted by column)
        [c.edge_c_csorted, csort_idx] = sort(c.edge_c);
        c.edge_r_csorted = c.edge_r(csort_idx);
        c.csort_idx = csort_idx;
        c.inv_csort_idx = zeros(c.ne, 1);
        c.inv_csort_idx(csort_idx) = (1:c.ne)';
        
        c.col_start = zeros(c.N, 1);
        c.col_end = zeros(c.N, 1);
        cur = 1;
        for ni = 1:c.N
            c.col_start(ni) = cur;
            while cur <= c.ne && c.edge_c_csorted(cur) == ni
                cur = cur + 1;
            end
            c.col_end(ni) = cur - 1;
        end
        
        % Find max row degree for padded vectorization
        row_deg = c.row_end - c.row_start + 1;
        c.max_row_deg = max(row_deg);
        
        cache = c;
    end
    
    % Retrieve from cache
    N = cache.N; M = cache.M;
    K_orig = cache.K_orig; K_ldpc = cache.K_ldpc;
    parity_len = cache.parity_len;
    ne = cache.ne; H = cache.H;
    edge_r = cache.edge_r; edge_c = cache.edge_c;

    % Setup LLRs — place info bits and parity bits in correct positions
    llr_full = zeros(N, 1);
    if K_orig > 0
        llr_full(1:K_orig) = llr_in(1:min(K_orig, length(llr_in)));
    end
    if parity_len > 0 && K_orig < length(llr_in)
        pl = min(parity_len, length(llr_in) - K_orig);
        llr_full(K_ldpc+1:K_ldpc+pl) = llr_in(K_orig+1:K_orig+pl);
    end
    % Filler bits are KNOWN to be zero → set large positive LLR
    % (positive LLR = bit is 0 with high confidence)
    num_filler = K_ldpc - K_orig;
    if num_filler > 0
        llr_full(K_orig+1:K_ldpc) = 1000;
    end

    % Initialize: r2c = 0, c2r = channel LLR for each edge
    r2c = zeros(ne, 1);
    c2r = llr_full(edge_c);

    success = false;
    alpha = 0.75;

    row_start = cache.row_start;
    row_end = cache.row_end;
    col_start = cache.col_start;
    col_end = cache.col_end;
    sort_idx = cache.sort_idx;
    inv_sort_idx = cache.inv_sort_idx;
    csort_idx = cache.csort_idx;
    inv_csort_idx = cache.inv_csort_idx;
    edge_c_sorted = cache.edge_c_sorted;

    for iter = 1:max_iter
        % ===== CHECK NODE UPDATE (Min-Sum) =====
        % Process each row: find min1, min2, product of signs
        r2c_new = zeros(ne, 1);
        for mi = 1:M
            rs = row_start(mi); re = row_end(mi);
            if rs > re, continue; end
            
            eidx = sort_idx(rs:re);  % original edge indices for this row
            vals = c2r(eidx);
            
            av = abs(vals);
            sg = sign(vals);
            sg(sg == 0) = 1;
            ps = prod(sg);
            
            % Find min and second-min
            [m1, mi1] = min(av);
            av2 = av; av2(mi1) = inf;
            m2 = min(av2);
            
            % Vectorized output computation
            mags = repmat(m1, length(eidx), 1);
            mags(mi1) = m2;
            r2c_new(eidx) = alpha * ps * sg .* mags;
        end
        r2c = r2c_new;

        % ===== VARIABLE NODE UPDATE =====
        hd = zeros(N, 1);
        for ni = 1:N
            cs = col_start(ni); ce = col_end(ni);
            if cs > ce
                hd(ni) = (llr_full(ni) < 0);
                continue;
            end
            eidx = csort_idx(cs:ce);
            tl = llr_full(ni) + sum(r2c(eidx));
            hd(ni) = (tl < 0);
            c2r(eidx) = tl - r2c(eidx);
        end

        % Syndrome check (every iteration for early termination)
        if mod(iter, 2) == 0 || iter == max_iter
            syn = mod(H * hd, 2);
            if all(syn == 0)
                success = true; break;
            end
        end
    end
    decoded = hd(1:K_orig);
end
