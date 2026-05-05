function coded = ldpc_encode_block(info_bits)
    persistent enc_cache
    
    info_bits = info_bits(:);
    K = length(info_bits);
    if K > 3824
        Kb = 22; nb = 68; mb = 46; bg = 1;
    else
        Kb = 10; nb = 52; mb = 42; bg = 2;
    end
    Z_set = sort([2 4 8 16 32 64 128 256 3 6 12 24 48 96 192 384 ...
                   5 10 20 40 80 160 320 7 14 28 56 112 224 ...
                   9 18 36 72 144 288 11 22 44 88 176 352 ...
                   13 26 52 104 208 15 30 60 120 240]);
    idx = find(Kb * Z_set >= K, 1);
    if isempty(idx), Z = Z_set(end); else, Z = Z_set(idx); end
    K_ldpc = Kb * Z;
    N = nb * Z;
    
    % Cache: build H once and precompute inv(Hp) for fast encoding
    if isempty(enc_cache) || enc_cache.K ~= K
        enc_cache.K = K;
        H = build_H(nb, mb, Z, bg);
        enc_cache.K_ldpc = K_ldpc;
        enc_cache.N = N;
        enc_cache.Hi = H(:, 1:K_ldpc);
        Hp = H(:, K_ldpc+1:end);
        
        % Compute inv(Hp) over GF(2) via Gaussian elimination on [Hp | I]
        [m, n] = size(Hp);
        fprintf('  LDPC init: computing inv(Hp) [%d x %d]... ', m, n);
        tic;
        Ab = [full(Hp) eye(m)];  % Augment with identity
        Ab = uint8(mod(Ab, 2));  % Use uint8 for fast bitxor
        pc = zeros(1, m);
        col = 1;
        for row = 1:m
            if col > n, break; end
            pr = 0;
            for r = row:m
                if Ab(r, col) == 1, pr = r; break; end
            end
            if pr == 0, col = col + 1; continue; end
            if pr ~= row
                Ab([row pr], :) = Ab([pr row], :);
            end
            pc(row) = col;
            % Vectorized row elimination using bitxor (much faster than mod)
            mask = (Ab(:, col) == 1);
            mask(row) = false;
            Ab(mask, :) = bitxor(Ab(mask, :), repmat(Ab(row, :), sum(mask), 1));
            col = col + 1;
        end
        % Extract inv(Hp) from the right side of the augmented matrix
        Hp_inv = double(Ab(:, n+1:end));
        enc_cache.Hp_inv = sparse(mod(Hp_inv, 2));
        fprintf('done (%.1fs)\n', toc);
    end
    
    num_filler = K_ldpc - K;
    info_padded = [info_bits; zeros(num_filler, 1)];
    
    % Fast encoding: s = Hi*info, parity = inv(Hp)*s
    s = mod(enc_cache.Hi * info_padded, 2);
    parity = mod(enc_cache.Hp_inv * s, 2);
    
    codeword = [info_padded; zeros(N - K_ldpc, 1)];
    codeword(K_ldpc+1:end) = parity;
    coded = [codeword(1:K); codeword(K_ldpc+1:end)];
end

function H = build_H(nb, mb, Z, bg)
    BM = -ones(mb, nb);
    if bg == 2
        entries = bg2_entries(Z);
        Kb = 10;
    else
        entries = bg1_entries(Z);
        Kb = 22;
    end
    for k = 1:size(entries, 1)
        BM(entries(k,1)+1, entries(k,2)+1) = mod(entries(k,3), Z);
    end
    % Generate extension rows (rows 4..mb-1 in 0-indexed, i.e. rows 5..mb in 1-indexed)
    for r = 5:mb  % 1-indexed row
        % Two info column connections with deterministic cyclic shifts
        sc1 = mod(r-5, Kb) + 1;           % 1-indexed info column
        BM(r, sc1) = mod(r*3, Z);
        sc2 = mod(r+3, Kb) + 1;           % 1-indexed info column
        BM(r, sc2) = mod(r*7+1, Z);
        % Parity column (identity circulant)
        pc = Kb + r;                       % 1-indexed parity column
        if pc <= nb, BM(r, pc) = 0; end
        % Staircase: connect to previous parity column
        if pc-1 > Kb && pc-1 <= nb, BM(r, pc-1) = mod(r, Z); end
    end
    N = nb * Z; M = mb * Z;
    H = sparse(M, N);
    for r = 1:mb
        for c = 1:nb
            sh = BM(r, c);
            if sh >= 0
                for i = 0:Z-1
                    j = mod(i + sh, Z);
                    H((r-1)*Z+i+1, (c-1)*Z+j+1) = 1;
                end
            end
        end
    end
end

function e = bg2_entries(~)
    % Base Graph 2 core rows (high connectivity rows 0-3)
    % Format: [row col shift]  (0-indexed, shift values are representative)
    e = [
        0 0 0; 0 1 0; 0 2 0; 0 3 0; 0 6 0; 0 9 0; 0 10 0; 0 11 0;
        1 0 0; 1 2 0; 1 3 0; 1 4 0; 1 7 0; 1 10 0; 1 12 0;
        2 0 0; 2 1 0; 2 4 0; 2 5 0; 2 8 0; 2 11 0; 2 13 0;
        3 1 0; 3 3 0; 3 5 0; 3 6 0; 3 9 0; 3 12 0; 3 14 0;
    ];
end

function e = bg1_entries(~)
    % Base Graph 1 core rows (high connectivity rows 0-3)
    % Format: [row col shift]  (0-indexed)
    e = [
        0 0 0; 0 1 0; 0 2 0; 0 3 0; 0 5 0; 0 7 0; 0 10 0; 0 12 0; 0 15 0; 0 18 0; 0 22 0; 0 23 0;
        1 0 0; 1 2 0; 1 4 0; 1 6 0; 1 8 0; 1 11 0; 1 13 0; 1 16 0; 1 19 0; 1 23 0; 1 24 0;
        2 1 0; 2 3 0; 2 5 0; 2 7 0; 2 9 0; 2 12 0; 2 14 0; 2 17 0; 2 20 0; 2 24 0; 2 25 0;
        3 0 0; 3 4 0; 3 6 0; 3 8 0; 3 11 0; 3 13 0; 3 15 0; 3 18 0; 3 21 0; 3 25 0; 3 26 0;
    ];
end
