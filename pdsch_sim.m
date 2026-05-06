function pdsch_sim()
%% 5G NR PDSCH Simulation: LDPC vs Convolutional Codes (K=7, K=9)
%  OFDM over L-tap Rayleigh fading channel, ZF & MMSE equalization
%  Project: EE5770 Codes and Waveforms

%% ============ BASE CONFIGURATION ============
cfg.A       = 1000;
cfg.RNTI    = 1;
cfg.cellId  = 0;
cfg.R       = 1/2;
cfg.M       = 2; % Bits per symbol: 2 (QPSK), 4 (16-QAM), 6 (64-QAM)
cfg.NFFT    = 256;
cfg.L       = 4; % number of channel taps
cfg.NCP     = cfg.L;
cfg.EbN0_dB = 0:1:8;
cfg.num_frames = 500;
cfg.ldpc_max_iter = 8;
cfg.conv_K  = 7;
cfg.conv_gen = [133 171];  % Octal: (133,171)_8
cfg.coding  = 'conv';
cfg.equalizer = 'mmse';
cfg.use_global_snr = false;  % false=per-subcarrier LLR, true=global SNR

% E is computed inside run_sim per coding scheme

%% ============ EXPERIMENT 1: Core BER/BLER (all schemes, L=4) ============
disp('========== EXPERIMENT 1: Core BER/BLER ==========');

disp('--- Conv K=7 / MMSE ---');
[ber_ck7m, bler_ck7m, dt_ck7m] = run_sim(cfg);

disp('--- Conv K=7 / ZF ---');
c2 = cfg; c2.equalizer = 'zf';
[ber_ck7z, bler_ck7z, dt_ck7z] = run_sim(c2);

disp('--- Conv K=9 / MMSE ---');
c3 = cfg; c3.conv_K = 9; c3.conv_gen = [561 753];
[ber_ck9m, bler_ck9m, dt_ck9m] = run_sim(c3);

disp('--- Conv K=9 / ZF ---');
c4 = cfg; c4.conv_K = 9; c4.conv_gen = [561 753]; c4.equalizer = 'zf';
[ber_ck9z, bler_ck9z, dt_ck9z] = run_sim(c4);

disp('--- LDPC / MMSE ---');
c5 = cfg; c5.coding = 'ldpc';
[ber_lm, bler_lm, dt_lm] = run_sim(c5);

disp('--- LDPC / ZF ---');
c6 = cfg; c6.coding = 'ldpc'; c6.equalizer = 'zf';
[ber_lz, bler_lz, dt_lz] = run_sim(c6);

% Plot 1: BER & BLER (all schemes)
figure('Name','Fig1: BER and BLER - All Schemes');
subplot(2,1,1);
semilogy(cfg.EbN0_dB, ber_ck7m, 'bo-', cfg.EbN0_dB, ber_ck7z, 'b^--', ...
         cfg.EbN0_dB, ber_ck9m, 'go-', cfg.EbN0_dB, ber_ck9z, 'g^--', ...
         cfg.EbN0_dB, ber_lm,   'rs-', cfg.EbN0_dB, ber_lz,   'r^--', 'LineWidth',1.5);
grid on; xlabel('E_b/N_0 (dB)'); ylabel('BER');
legend('Conv K=7 MMSE','Conv K=7 ZF','Conv K=9 MMSE','Conv K=9 ZF', ...
       'LDPC MMSE','LDPC ZF','Location','southwest');
title('BER: LDPC vs Conv K=7 vs Conv K=9');
subplot(2,1,2);
semilogy(cfg.EbN0_dB, max(bler_ck7m,1e-4), 'bo-', cfg.EbN0_dB, max(bler_ck7z,1e-4), 'b^--', ...
         cfg.EbN0_dB, max(bler_ck9m,1e-4), 'go-', cfg.EbN0_dB, max(bler_ck9z,1e-4), 'g^--', ...
         cfg.EbN0_dB, max(bler_lm,1e-4),   'rs-', cfg.EbN0_dB, max(bler_lz,1e-4),   'r^--', 'LineWidth',1.5);
grid on; xlabel('E_b/N_0 (dB)'); ylabel('BLER');
legend('Conv K=7 MMSE','Conv K=7 ZF','Conv K=9 MMSE','Conv K=9 ZF', ...
       'LDPC MMSE','LDPC ZF','Location','southwest');
title('BLER: LDPC vs Conv K=7 vs Conv K=9');

%% ============ EXPERIMENT 2: ZF vs MMSE at L=4 and L=16 ============
disp('========== EXPERIMENT 2: ZF vs MMSE, L=4 and L=16 ==========');

disp('--- LDPC MMSE L=16 ---');
c7 = cfg; c7.coding = 'ldpc'; c7.L = 16; c7.NCP = 16;
[ber_lm16, bler_lm16, ~] = run_sim(c7);
disp('--- LDPC ZF L=16 ---');
c8 = cfg; c8.coding = 'ldpc'; c8.L = 16; c8.NCP = 16; c8.equalizer = 'zf';
[ber_lz16, bler_lz16, ~] = run_sim(c8);

figure('Name','Fig2: ZF vs MMSE (L=4 and L=16)');
subplot(2,1,1);
semilogy(cfg.EbN0_dB, ber_lm,   'bs-',  cfg.EbN0_dB, ber_lz,   'b^--', ...
         cfg.EbN0_dB, ber_lm16, 'rs-',  cfg.EbN0_dB, ber_lz16, 'r^--', 'LineWidth',1.5);
grid on; xlabel('E_b/N_0 (dB)'); ylabel('BER');
legend('MMSE L=4','ZF L=4','MMSE L=16','ZF L=16','Location','southwest');
title('BER: ZF vs MMSE Equalizer');
subplot(2,1,2);
semilogy(cfg.EbN0_dB, max(bler_lm,1e-4),   'bs-',  cfg.EbN0_dB, max(bler_lz,1e-4),   'b^--', ...
         cfg.EbN0_dB, max(bler_lm16,1e-4), 'rs-',  cfg.EbN0_dB, max(bler_lz16,1e-4), 'r^--', 'LineWidth',1.5);
grid on; xlabel('E_b/N_0 (dB)'); ylabel('BLER');
legend('MMSE L=4','ZF L=4','MMSE L=16','ZF L=16','Location','southwest');
title('BLER: ZF vs MMSE Equalizer');

%% ============ EXPERIMENT 3: Per-subcarrier SNR CDF ============
disp('========== EXPERIMENT 3: SNR Distribution CDF ==========');
[snr_zf, snr_mmse] = collect_snr_cdf(cfg, 5);
figure('Name','Fig3: Per-Subcarrier SNR CDF');
cdfplot(10*log10(max(snr_zf, 1e-10))); hold on;
cdfplot(10*log10(max(snr_mmse, 1e-10)));
xlabel('\gamma_k (dB)'); ylabel('CDF'); grid on;
legend('ZF','MMSE','Location','southeast');
title('Per-Subcarrier SNR Distribution (E_b/N_0 = 5 dB, L=4)');

%% ============ EXPERIMENT 4: LLR Mismatch ============
disp('========== EXPERIMENT 4: LLR Mismatch Experiment ==========');
disp('--- LDPC MMSE (global SNR) ---');
c9 = cfg; c9.coding = 'ldpc'; c9.use_global_snr = true;
[~, bler_lm_g, ~] = run_sim(c9);
disp('--- LDPC ZF (global SNR) ---');
c10 = cfg; c10.coding = 'ldpc'; c10.equalizer = 'zf'; c10.use_global_snr = true;
[~, bler_lz_g, ~] = run_sim(c10);

figure('Name','Fig4: LLR Mismatch');
semilogy(cfg.EbN0_dB, max(bler_lm,1e-4),  'bs-',  cfg.EbN0_dB, max(bler_lz,1e-4),  'b^--', ...
         cfg.EbN0_dB, max(bler_lm_g,1e-4), 'rs-', cfg.EbN0_dB, max(bler_lz_g,1e-4), 'r^--', 'LineWidth',1.5);
grid on; xlabel('E_b/N_0 (dB)'); ylabel('BLER');
legend('MMSE (per-SC)','ZF (per-SC)','MMSE (global)','ZF (global)','Location','southwest');
title('LLR Mismatch: Per-Subcarrier vs Global SNR Scaling');

%% ============ EXPERIMENT 5: Decoding Complexity ============
figure('Name','Fig5: Decoding Complexity');
subplot(1,2,1);
bar(categorical({'Conv K=7','Conv K=9','LDPC'}), [mean(dt_ck7m) mean(dt_ck9m) mean(dt_lm)]);
ylabel('Avg Decode Time (s)'); title('Mean Decode Time (MMSE)');
subplot(1,2,2);
plot(cfg.EbN0_dB, dt_ck7m, 'bo-', cfg.EbN0_dB, dt_ck9m, 'go-', cfg.EbN0_dB, dt_lm, 'rs-', 'LineWidth',1.5);
xlabel('E_b/N_0 (dB)'); ylabel('Avg Decode Time (s)');
legend('Conv K=7','Conv K=9','LDPC'); grid on; title('Decode Time vs SNR');

disp('========== ALL EXPERIMENTS COMPLETE ==========');
end

%% ================================================================
function [ber_res, bler_res, dt_res] = run_sim(cfg)
    ber_res  = zeros(size(cfg.EbN0_dB));
    bler_res = zeros(size(cfg.EbN0_dB));
    dt_res   = zeros(size(cfg.EbN0_dB));
    % Compute E based on coding scheme
    B = cfg.A + 24;  % Info + CRC24A
    if strcmpi(cfg.coding, 'ldpc')
        % LDPC: E determined by target rate (LDPC has internal rate matching)
        cfg.E = cfg.M * ceil(cfg.A / cfg.R / cfg.M);
    else
        % Conv: E = full coded length (no puncturing, tail included)
        coded_len_conv = length(cfg.conv_gen) * (B + cfg.conv_K - 1);
        cfg.E = cfg.M * ceil(coded_len_conv / cfg.M);
    end
    Reff = B / cfg.E;

    for si = 1:length(cfg.EbN0_dB)
        EbN0 = cfg.EbN0_dB(si);
        sigma2_w = 1 / (2 * Reff * cfg.M * 10^(EbN0/10));
        tbe = 0; tb = 0; tble = 0; tbl = 0; tdt = 0;

        for fr = 1:cfg.num_frames
            % ---------------- TRANSMITTER ----------------
            tx_bits = randi([0 1], cfg.A, 1);
            tx_crc = crc24a_encode(tx_bits);
            code_blocks = code_block_segment(tx_crc, cfg.coding);
            
            coded_bits = [];
            cb_coded_lens = zeros(length(code_blocks),1);
            for cb = 1:length(code_blocks)
                if strcmpi(cfg.coding, 'ldpc')
                    cc = ldpc_encode_block(code_blocks{cb});
                else
                    cc = conv_encode(code_blocks{cb}, cfg.conv_K, cfg.conv_gen);
                end
                cb_coded_lens(cb) = length(cc);
                coded_bits = [coded_bits; cc];
            end

            rm_bits = rate_match(coded_bits, cfg.E);
            scr_seq = scramble_seq(cfg.E, cfg.RNTI, cfg.cellId);
            scr_bits = xor(rm_bits(:), scr_seq(:));
            symbols = gray_map(scr_bits, cfg.M);
            
            [tx_sig, sc_idx, num_sym] = ofdm_modulate(symbols, cfg.NFFT, cfg.NCP);
            
            % ---------------- CHANNEL ----------------
            [rx_sig, h_taps] = channel_multipath(tx_sig, cfg.L, cfg.NFFT, cfg.NCP, num_sym, sigma2_w);
            
            % ---------------- RECEIVER ----------------
            [eq_sym, snr_sc] = ofdm_demod_equalize(rx_sig, h_taps, cfg.NFFT, cfg.NCP, num_sym, sc_idx, sigma2_w, cfg.equalizer);

            eq_sym = eq_sym(1:length(symbols));
            snr_sc = snr_sc(1:length(symbols));

            % LLR mismatch experiment: use global (mean) SNR if requested
            if isfield(cfg, 'use_global_snr') && cfg.use_global_snr
                snr_sc(:) = mean(snr_sc);
            end

            llr = gray_demap_llr(eq_sym, cfg.M, snr_sc);
            llr_d = descramble_llr(llr, scr_seq);
            
            % RX side calculation of expected coded length (NO CHEATING!)
            B_rx = cfg.A + 24; 
            if strcmpi(cfg.coding, 'ldpc')
                Kcb_rx = 8448;
            else
                Kcb_rx = 6144;
            end
            C_rx = ceil(B_rx / (Kcb_rx - 24));
            seg_len_rx = ceil(B_rx / C_rx);
            
            expected_coded_len = 0;
            cb_lens_rx = zeros(C_rx,1);
            il_rx = zeros(C_rx,1); 
            
            for c_rx = 1:C_rx
                si_rx = (c_rx-1)*seg_len_rx + 1;
                ei_rx = min(c_rx*seg_len_rx, B_rx);
                
                % CONDITIONAL CRC24B: Only attach if multiple segments exist
                if C_rx > 1
                    seg_bits = (ei_rx - si_rx + 1) + 24; 
                else
                    seg_bits = (ei_rx - si_rx + 1);      
                end
                il_rx(c_rx) = seg_bits; 
                
                if strcmpi(cfg.coding, 'ldpc')
                    if seg_bits > 3824, Kb_rx = 22; nb_rx = 68; else, Kb_rx = 10; nb_rx = 52; end
                    Z_set = sort([2 4 8 16 32 64 128 256 3 6 12 24 48 96 192 384 5 10 20 40 80 160 320 7 14 28 56 112 224 9 18 36 72 144 288 11 22 44 88 176 352 13 26 52 104 208 15 30 60 120 240]);
                    idx_Z = find(Kb_rx * Z_set >= seg_bits, 1);
                    if isempty(idx_Z), Z_rx = Z_set(end); else, Z_rx = Z_set(idx_Z); end
                    K_ldpc_rx = Kb_rx * Z_rx;
                    N_rx = nb_rx * Z_rx;
                    num_filler_rx = K_ldpc_rx - seg_bits;
                    % Coded output = info(K) + parity(N-K_ldpc), filler excluded
                    cb_lens_rx(c_rx) = seg_bits + (N_rx - K_ldpc_rx);
                    expected_coded_len = expected_coded_len + cb_lens_rx(c_rx);
                else
                    cb_lens_rx(c_rx) = 2 * (seg_bits + cfg.conv_K - 1);
                    expected_coded_len = expected_coded_len + cb_lens_rx(c_rx);
                end
            end
            
            llr_rec = rate_recover(llr_d, expected_coded_len, cfg.E);

            % Timer starts exactly here for accurate matched-complexity
            tic;
            dec_blocks = {};
            blk_pass = true;
            ptr = 1;
            
            for cb = 1:C_rx
                cl = cb_lens_rx(cb);
                cb_llr = llr_rec(ptr:ptr+cl-1);
                ptr = ptr + cl;
                
                il = il_rx(cb); 
                
                if strcmpi(cfg.coding, 'ldpc')
                    [db, ~] = ldpc_decode_block(cb_llr, cfg.ldpc_max_iter, il);
                else
                    db = viterbi_decode(cb_llr, cfg.conv_K, cfg.conv_gen);
                    db = db(1:min(il,length(db)));
                    if length(db) < il
                        db = [db; zeros(il-length(db),1)];
                    end
                end
                
                if C_rx > 1
                    [dp, cok] = crc24b_check(db);
                    if ~cok, blk_pass = false; end
                    dec_blocks{cb} = dp;
                else
                    dec_blocks{cb} = db;
                end
            end
            tdt = tdt + toc;

            rx_payload = [];
            for cb = 1:length(dec_blocks)
                rx_payload = [rx_payload; dec_blocks{cb}(:)];
            end
            [rx_data, crc_pass] = crc24a_check(rx_payload);
            rx_data = rx_data(1:min(length(rx_data), cfg.A));
            clen = min(length(rx_data), cfg.A);
            if clen < cfg.A
                rx_data = [rx_data; zeros(cfg.A-clen,1)];
            end
            tbe = tbe + sum(rx_data ~= tx_bits);
            tb = tb + cfg.A;
            if ~crc_pass || ~blk_pass
                tble = tble + 1;
            end
            tbl = tbl + 1;
        end
        ber_res(si) = tbe/tb;
        bler_res(si) = tble/tbl;
        dt_res(si) = tdt/cfg.num_frames;
        fprintf('%s/%s Eb/N0=%2d dB: BER=%.2e BLER=%.4f DecTime=%.4fs\n', ...
            cfg.coding, cfg.equalizer, EbN0, ber_res(si), bler_res(si), dt_res(si));
    end
end

%% ==================== SNR CDF COLLECTION ====================
function [snr_zf_all, snr_mmse_all] = collect_snr_cdf(cfg, EbN0_fixed)
    % Collect per-subcarrier SNR samples for CDF plot at fixed Eb/N0
    B = cfg.A + 24;
    if strcmpi(cfg.coding, 'ldpc')
        E_tmp = cfg.M * ceil(cfg.A / cfg.R / cfg.M);
    else
        E_tmp = cfg.M * ceil(length(cfg.conv_gen) * (B + cfg.conv_K - 1) / cfg.M);
    end
    Reff = B / E_tmp;
    sigma2_w = 1 / (2 * Reff * cfg.M * 10^(EbN0_fixed/10));
    
    num_data_sc = floor(cfg.NFFT * 0.75);
    num_data_sc = num_data_sc - mod(num_data_sc, 2);
    sc_idx = [2:num_data_sc/2+1, cfg.NFFT-num_data_sc/2+1:cfg.NFFT].';
    
    n_frames = min(cfg.num_frames, 200);
    snr_zf_all = [];
    snr_mmse_all = [];
    
    for fr = 1:n_frames
        % Generate channel taps
        h = (randn(cfg.L, 1) + 1j*randn(cfg.L, 1)) / sqrt(2*cfg.L);
        hp = [h; zeros(cfg.NFFT - cfg.L, 1)];
        H = fft(hp);
        Hd = H(sc_idx);
        
        % ZF SNR
        gam_zf = abs(Hd).^2 / sigma2_w;
        % MMSE SNR
        den = abs(Hd).^2 + sigma2_w;
        gam_mmse = den ./ sigma2_w;
        
        snr_zf_all = [snr_zf_all; gam_zf];
        snr_mmse_all = [snr_mmse_all; gam_mmse];
    end
    fprintf('  Collected %d SNR samples for CDF\n', length(snr_zf_all));
end

%% ==================== CRC ====================
function crc = crc_calc(bits, poly)
    bits = bits(:).';
    crc_len = length(poly) - 1;
    reg = zeros(1, crc_len);
    for i = 1:length(bits)
        fb = xor(bits(i), reg(1));
        reg = [reg(2:end) 0];
        if fb
            reg = xor(reg, poly(2:end));
        end
    end
    crc = reg(:);
end

function out = crc24a_encode(bits)
    poly = [1 1 0 0 0 0 1 1 0 0 1 0 0 1 1 0 0 1 1 1 1 1 0 1 1];
    out = [bits(:); crc_calc(bits, poly)];
end

function [payload, pass] = crc24a_check(bits)
    bits = bits(:);
    if length(bits) < 25, payload = bits; pass = false; return; end
    payload = bits(1:end-24);
    rcrc = bits(end-23:end);
    poly = [1 1 0 0 0 0 1 1 0 0 1 0 0 1 1 0 0 1 1 1 1 1 0 1 1];
    pass = isequal(rcrc(:), crc_calc(payload, poly));
end

function out = crc24b_encode(bits)
    poly = [1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 1 1 1];
    out = [bits(:); crc_calc(bits, poly)];
end

function [payload, pass] = crc24b_check(bits)
    bits = bits(:);
    if length(bits) < 25, payload = bits; pass = false; return; end
    payload = bits(1:end-24);
    rcrc = bits(end-23:end);
    poly = [1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 1 1 1];
    pass = isequal(rcrc(:), crc_calc(payload, poly));
end

%% ==================== CODE BLOCK SEGMENTATION ====================
function blocks = code_block_segment(bits, coding_type)
    bits = bits(:);
    B = length(bits);
    if strcmpi(coding_type, 'ldpc')
        Kcb = 8448;
    else
        Kcb = 6144;
    end
    if B <= Kcb
        blocks = {bits};
    else
        C = ceil(B / (Kcb - 24));
        seg_len = ceil(B / C);
        blocks = cell(1, C);
        for c = 1:C
            si = (c-1)*seg_len + 1;
            ei = min(c*seg_len, B);
            blocks{c} = crc24b_encode(bits(si:ei));
        end
    end
end

%% ==================== LDPC ENCODER ====================
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
        0 0 307; 0 1 19; 0 2 50; 0 3 369; 0 5 181; 0 6 216; 0 10 288; 0 12 17; 0 15 215; 0 18 242; 0 22 1; 0 23 0;
        1 0 76; 1 2 76; 1 4 288; 1 7 331; 1 8 331; 1 11 295; 1 12 342; 1 16 354; 1 19 331; 1 23 0; 1 24 0;
        2 1 250; 2 4 332; 2 5 256; 2 7 267; 2 9 63; 2 10 129; 2 13 200; 2 17 131; 2 20 13; 2 24 0; 2 25 0;
        3 0 276; 3 4 275; 3 6 199; 3 8 56; 3 11 305; 3 13 341; 3 17 300; 3 18 271; 3 21 357; 3 22 1; 3 25 0;
    ];
end

%% ==================== VECTORIZED LDPC DECODER ====================
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

%% ==================== CONVOLUTIONAL ENCODER ====================
function coded = conv_encode(info_bits, K, gen_poly)
    info_bits = info_bits(:).';
    ng = length(gen_poly);
    mem = K - 1;
    gb = zeros(ng, K);
    for g = 1:ng
        % gen_poly values like 133, 171 are OCTAL representations stored as decimal
        % Use num2str to get the digit string directly (each digit IS an octal digit)
        bs = num2str(gen_poly(g));
        bv = [];
        for ch = 1:length(bs)
            d = str2double(bs(ch));
            bv = [bv bitget(d, 3:-1:1)];
        end
        if length(bv) >= K
            gb(g,:) = bv(end-K+1:end);
        else
            gb(g, K-length(bv)+1:end) = bv;
        end
    end
    iwt = [info_bits zeros(1, mem)];
    L = length(iwt);
    coded = zeros(ng * L, 1);
    state = zeros(1, mem);
    for i = 1:L
        reg = [iwt(i) state];
        for g = 1:ng
            coded(ng*(i-1)+g) = mod(sum(reg .* gb(g,:)), 2);
        end
        state = [iwt(i) state(1:end-1)];
    end
end

%% ==================== VITERBI DECODER ====================
function decoded = viterbi_decode(llr, K, gen_poly)
    persistent cached_K cached_gen_poly cached_nxt cached_outp cached_ns cached_ng cached_mem
    
    llr = llr(:).';
    L = floor(length(llr) / length(gen_poly));
    
    if isempty(cached_K) || cached_K ~= K || ~isequal(cached_gen_poly, gen_poly)
        cached_K = K; 
        cached_gen_poly = gen_poly;
        cached_ng = length(gen_poly);
        cached_mem = K - 1;
        cached_ns = 2^cached_mem;
        
        gb = zeros(cached_ng, K);
        for g = 1:cached_ng
            % gen_poly values like 133, 171 are OCTAL representations stored as decimal
            bs = num2str(gen_poly(g));
            bv = [];
            for ch = 1:length(bs)
                d = str2double(bs(ch));
                bv = [bv bitget(d, 3:-1:1)];
            end
            if length(bv) >= K
                gb(g,:) = bv(end-K+1:end);
            else
                gb(g, K-length(bv)+1:end) = bv;
            end
        end
        
        cached_nxt = zeros(cached_ns, 2);
        cached_outp = zeros(cached_ns, 2, cached_ng);
        for s = 0:cached_ns-1
            sb = bitget(s, cached_mem:-1:1);
            for inp = 0:1
                reg = [inp sb];
                nsb = reg(1:cached_mem);
                cached_nxt(s+1, inp+1) = 0;
                for b = 1:cached_mem
                    cached_nxt(s+1, inp+1) = cached_nxt(s+1, inp+1) + nsb(b) * 2^(cached_mem-b);
                end
                for g = 1:cached_ng
                    cached_outp(s+1, inp+1, g) = mod(sum(reg .* gb(g,:)), 2);
                end
            end
        end
    end
    
    nxt = cached_nxt;
    outp = cached_outp;
    ns = cached_ns;
    ng = cached_ng;
    mem = cached_mem;

    pm = -inf(ns, 1);
    pm(1) = 0;
    surv = zeros(ns, L);
    sh = zeros(ns, L);
    
    for t = 1:L
        npm = -inf(ns, 1);
        nsu = zeros(ns, 1);
        nsh = zeros(ns, 1);
        bl = llr((t-1)*ng+1 : t*ng);
        for s = 0:ns-1
            if pm(s+1) == -inf, continue; end
            for inp = 0:1
                nst = nxt(s+1, inp+1);
                bm = 0;
                for g = 1:ng
                    ex = 1 - 2*outp(s+1, inp+1, g);
                    bm = bm + bl(g) * ex;
                end
                cand = pm(s+1) + bm;
                if cand > npm(nst+1)
                    npm(nst+1) = cand;
                    nsu(nst+1) = inp;
                    nsh(nst+1) = s;
                end
            end
        end
        pm = npm;
        surv(:,t) = nsu;
        sh(:,t) = nsh;
    end
    decoded = zeros(L, 1);
    s = 0;
    for t = L:-1:1
        decoded(t) = surv(s+1, t);
        s = sh(s+1, t);
    end
    decoded = decoded(1:L-mem);
end

%% ==================== RATE MATCHING ====================
function out = rate_match(coded, E)
    coded = coded(:); N = length(coded);
    if E <= N
        out = coded(1:E);
    else
        out = zeros(E, 1);
        for i = 1:E
            out(i) = coded(mod(i-1, N) + 1);
        end
    end
end

function llr_out = rate_recover(llr_in, coded_len, E)
    llr_in = llr_in(:);
    if E <= coded_len
        llr_out = zeros(coded_len, 1);
        llr_out(1:E) = llr_in;
    else
        llr_out = zeros(coded_len, 1);
        for i = 1:E
            idx = mod(i-1, coded_len) + 1;
            llr_out(idx) = llr_out(idx) + llr_in(i);
        end
    end
end

%% ==================== SCRAMBLING ====================
function seq = scramble_seq(len, rnti, cellId)
    Nc = 1600;
    c_init = rnti * 2^15 + cellId;
    tl = len + Nc;
    x1 = zeros(1, tl + 31);
    x2 = zeros(1, tl + 31);
    x1(1) = 1;
    ib = zeros(1, 31);
    tmp = c_init;
    for b = 1:31
        ib(b) = mod(tmp, 2);
        tmp = floor(tmp / 2);
    end
    x2(1:31) = ib;
    for n = 1:tl
        x1(n+31) = mod(x1(n+3) + x1(n), 2);
        x2(n+31) = mod(x2(n+3) + x2(n+2) + x2(n+1) + x2(n), 2);
    end
    seq = zeros(len, 1);
    for n = 1:len
        seq(n) = mod(x1(n+Nc) + x2(n+Nc), 2);
    end
end

function llr_out = descramble_llr(llr_in, scr_seq)
    llr_in = llr_in(:); scr_seq = scr_seq(:);
    sl = min(length(llr_in), length(scr_seq));
    llr_out = llr_in;
    llr_out(1:sl) = llr_in(1:sl) .* (1 - 2*scr_seq(1:sl));
end

%% ==================== MODULATION ====================
function symbols = gray_map(bits, M)
    bits = bits(:);
    ns = floor(length(bits) / M);
    bits = bits(1:ns*M);
    symbols = zeros(ns, 1);
    if M == 2
        for i = 1:ns
            b = bits((i-1)*2+1 : i*2);
            symbols(i) = (1-2*b(1))/sqrt(2) + 1j*(1-2*b(2))/sqrt(2);
        end
    elseif M == 4
        nm = 1/sqrt(10);
        for i = 1:ns
            b = bits((i-1)*4+1 : i*4);
            re = (1-2*b(1))*(2-(1-2*b(2)));
            im = (1-2*b(3))*(2-(1-2*b(4)));
            symbols(i) = nm*(re+1j*im);
        end
    elseif M == 6
        nm = 1/sqrt(42);
        for i = 1:ns
            b = bits((i-1)*6+1 : i*6);
            re = (1-2*b(1))*(4-(1-2*b(2))*(2-(1-2*b(3))));
            im = (1-2*b(4))*(4-(1-2*b(5))*(2-(1-2*b(6))));
            symbols(i) = nm*(re+1j*im);
        end
    end
end

function llr = gray_demap_llr(symbols, M, snr)
    symbols = symbols(:); snr = snr(:);
    ns = length(symbols);
    llr = zeros(ns*M, 1);
    if M == 2
        for i = 1:ns
            g = snr(i);
            llr((i-1)*2+1) = 4*g*real(symbols(i))/sqrt(2);
            llr((i-1)*2+2) = 4*g*imag(symbols(i))/sqrt(2);
        end
    else
        [con, bt] = get_constellation(M);
        for i = 1:ns
            g = snr(i);
            for b = 1:M
                md0 = inf; md1 = inf;
                for s = 1:length(con)
                    d = abs(symbols(i) - con(s))^2;
                    if bt(s,b) == 0
                        if d < md0, md0 = d; end
                    else
                        if d < md1, md1 = d; end
                    end
                end
                llr((i-1)*M+b) = g*(md1-md0);
            end
        end
    end
end

function [con, bt] = get_constellation(M)
    if M == 4
        nm = 1/sqrt(10);
        vals = [-3 -1 1 3];
        con = zeros(16,1); bt = zeros(16,4);
        idx = 1;
        for re = vals
            for im = vals
                con(idx) = nm*(re+1j*im);
                bt(idx,1) = (re<0);
                bt(idx,2) = (abs(re)>2);
                bt(idx,3) = (im<0);
                bt(idx,4) = (abs(im)>2);
                idx = idx+1;
            end
        end
    else
        nm = 1/sqrt(42);
        vals = [-7 -5 -3 -1 1 3 5 7];
        con = zeros(64,1); bt = zeros(64,6);
        idx = 1;
        for re = vals
            for im = vals
                con(idx) = nm*(re+1j*im);
                bt(idx,1) = (re<0);
                bt(idx,2) = (abs(re)>4);
                bt(idx,3) = (abs(abs(re)-4)>2);
                bt(idx,4) = (im<0);
                bt(idx,5) = (abs(im)>4);
                bt(idx,6) = (abs(abs(im)-4)>2);
                idx = idx+1;
            end
        end
    end
end

%% ==================== OFDM ====================
function [tx_sig, sc_idx, num_sym] = ofdm_modulate(symbols, NFFT, NCP)
    symbols = symbols(:);
    nds = floor(NFFT * 0.75);
    nds = nds - mod(nds, 2);
    num_sym = ceil(length(symbols) / nds);
    total_s = num_sym * nds;
    symbols = [symbols; zeros(total_s - length(symbols), 1)];
    half = nds / 2;
    sc_idx = [2:half+1, NFFT-half+1:NFFT].';
    
    sym_len = NFFT + NCP;
    tx_sig = zeros(num_sym * sym_len, 1);
    
    for s = 1:num_sym
        S = zeros(NFFT, 1);
        S(sc_idx) = symbols((s-1)*nds+1 : s*nds);
        x = ifft(S) * sqrt(NFFT);
        start_idx = (s-1)*sym_len + 1;
        tx_sig(start_idx : start_idx + sym_len - 1) = [x(end-NCP+1:end); x];
    end
end

function [eq_sym, snr_out] = ofdm_demod_equalize(rx_sig, h_taps, NFFT, NCP, num_sym, sc_idx, sigma2_w, eq_type)
    num_sc = length(sc_idx);
    eq_sym = zeros(num_sym * num_sc, 1);
    snr_out = zeros(num_sym * num_sc, 1);
    
    for s = 1:num_sym
        si = (s-1)*(NFFT+NCP) + 1;
        if si+NCP+NFFT-1 > length(rx_sig), break; end
        y = rx_sig(si+NCP : si+NCP+NFFT-1);
        Y = fft(y) / sqrt(NFFT);
        h = h_taps{s};
        hp = [h(:); zeros(NFFT-length(h), 1)];
        H = fft(hp);
        Yd = Y(sc_idx);
        Hd = H(sc_idx);
        
        if strcmpi(eq_type, 'zf')
            Sh = Yd ./ Hd;
            gam = abs(Hd).^2 / sigma2_w;
        else
            % Proper MMSE equalization: output = mu*X + n_eff
            den = abs(Hd).^2 + sigma2_w;
            Sh = conj(Hd) .* Yd ./ den;
            % Effective SNR for LLR: g = a/sigma2_n = (|H|^2+sigma2_w)/sigma2_w
            gam = den ./ sigma2_w;
        end
        
        idx_start = (s-1)*num_sc + 1;
        eq_sym(idx_start : idx_start + num_sc - 1) = Sh;
        snr_out(idx_start : idx_start + num_sc - 1) = gam;
    end
end

%% ==================== CHANNEL ====================
function [rx_sig, h_taps] = channel_multipath(tx_sig, L, NFFT, NCP, num_sym, sigma2_w)
    sl = NFFT + NCP;
    rx_sig = zeros(num_sym * sl, 1);
    h_taps = cell(num_sym, 1);
    
    for s = 1:num_sym
        si = (s-1)*sl + 1;
        ei = si + sl - 1;
        if ei > length(tx_sig)
            xs = [tx_sig(si:end); zeros(ei-length(tx_sig), 1)];
        else
            xs = tx_sig(si:ei);
        end
        h = (randn(L,1) + 1j*randn(L,1)) / sqrt(2*L);
        h_taps{s} = h;
        yc = conv(xs, h);
        ys = yc(1:sl);
        noise = sqrt(sigma2_w/2) * (randn(sl,1) + 1j*randn(sl,1));
        
        rx_sig((s-1)*sl+1 : s*sl) = ys + noise;
    end
end
