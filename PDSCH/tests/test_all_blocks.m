function test_all_blocks()
%% Unit Tests for PDSCH Simulation Processing Blocks
%  Run: >> test_all_blocks
%  All tests print PASS/FAIL. Errors halt execution.

fprintf('======== PDSCH Unit Tests ========\n\n');

test_crc();
test_code_block_segment();
test_conv_encode_decode();
test_ldpc_encode_decode();
test_ofdm();
test_equalizers();
test_rate_matching();
test_scrambling();

fprintf('\n======== ALL TESTS PASSED ========\n');
end

%% ==================== CRC Tests ====================
function test_crc()
    fprintf('--- CRC Tests ---\n');
    
    % CRC24A: encode then check must pass
    bits = randi([0 1], 100, 1);
    encoded = crc24a_encode(bits);
    assert(length(encoded) == 124, 'CRC24A length mismatch');
    [payload, pass] = crc24a_check(encoded);
    assert(pass, 'CRC24A check failed on valid data');
    assert(isequal(payload, bits), 'CRC24A payload mismatch');
    
    % Flip a bit — CRC must fail
    corrupted = encoded;
    corrupted(50) = 1 - corrupted(50);
    [~, pass2] = crc24a_check(corrupted);
    assert(~pass2, 'CRC24A should fail on corrupted data');
    
    % CRC24B: encode then check
    encoded_b = crc24b_encode(bits);
    assert(length(encoded_b) == 124, 'CRC24B length mismatch');
    [payload_b, pass_b] = crc24b_check(encoded_b);
    assert(pass_b, 'CRC24B check failed on valid data');
    assert(isequal(payload_b, bits), 'CRC24B payload mismatch');
    
    fprintf('  CRC24A encode/check: PASS\n');
    fprintf('  CRC24A corruption detect: PASS\n');
    fprintf('  CRC24B encode/check: PASS\n');
end

%% ==================== Code Block Segmentation ====================
function test_code_block_segment()
    fprintf('--- Code Block Segmentation Tests ---\n');
    
    % Small block — no segmentation
    bits = randi([0 1], 1000, 1);
    blocks = code_block_segment(bits, 'ldpc');
    assert(length(blocks) == 1, 'Should be 1 block for 1000 bits LDPC');
    assert(length(blocks{1}) == 1000, 'Block size mismatch');
    
    blocks_c = code_block_segment(bits, 'conv');
    assert(length(blocks_c) == 1, 'Should be 1 block for 1000 bits conv');
    
    % Large block — should segment
    big_bits = randi([0 1], 10000, 1);
    blocks_big = code_block_segment(big_bits, 'ldpc');
    assert(length(blocks_big) > 1, 'Should segment >8448 bits for LDPC');
    
    % Verify total bits preserved (each segment gets CRC24B)
    total = 0;
    for i = 1:length(blocks_big)
        total = total + length(blocks_big{i}) - 24; % subtract CRC24B
    end
    assert(total == 10000, 'Segmented bits total mismatch');
    
    fprintf('  Single block (LDPC): PASS\n');
    fprintf('  Single block (Conv): PASS\n');
    fprintf('  Multi-block segmentation: PASS\n');
end

%% ==================== Convolutional Encode/Decode ====================
function test_conv_encode_decode()
    fprintf('--- Convolutional Encoder/Decoder Tests ---\n');
    
    % K=7 test
    info = randi([0 1], 200, 1);
    K = 7; gen = [133 171];
    coded = conv_encode(info, K, gen);
    expected_len = 2 * (length(info) + K - 1);
    assert(length(coded) == expected_len, 'Conv K=7 coded length wrong');
    
    % Noiseless decode — must be perfect
    llr = 1 - 2*coded;  % Perfect LLRs (bit 0 -> +1, bit 1 -> -1)
    llr = llr * 100;     % Strong LLRs
    decoded = viterbi_decode(llr, K, gen);
    decoded = decoded(1:length(info));
    assert(isequal(decoded(:), info(:)), 'Conv K=7 noiseless decode failed');
    
    % K=9 test
    K9 = 9; gen9 = [561 753];
    coded9 = conv_encode(info, K9, gen9);
    expected_len9 = 2 * (length(info) + K9 - 1);
    assert(length(coded9) == expected_len9, 'Conv K=9 coded length wrong');
    
    llr9 = (1 - 2*coded9) * 100;
    decoded9 = viterbi_decode(llr9, K9, gen9);
    decoded9 = decoded9(1:length(info));
    assert(isequal(decoded9(:), info(:)), 'Conv K=9 noiseless decode failed');
    
    fprintf('  Conv K=7 encode length: PASS\n');
    fprintf('  Conv K=7 noiseless roundtrip: PASS\n');
    fprintf('  Conv K=9 encode length: PASS\n');
    fprintf('  Conv K=9 noiseless roundtrip: PASS\n');
end

%% ==================== LDPC Encode/Decode ====================
function test_ldpc_encode_decode()
    fprintf('--- LDPC Encoder/Decoder Tests ---\n');
    
    info = randi([0 1], 500, 1);
    coded = ldpc_encode_block(info);
    assert(length(coded) > length(info), 'LDPC coded should be longer');
    
    % Verify codeword satisfies H*c = 0 mod 2
    % (We test via noiseless decode)
    K = length(info);
    if K > 3824, Kb=22; nb=68; mb=46; else, Kb=10; nb=52; mb=42; end
    Z_set = sort([2 4 8 16 32 64 128 256 3 6 12 24 48 96 192 384 ...
                   5 10 20 40 80 160 320 7 14 28 56 112 224 ...
                   9 18 36 72 144 288 11 22 44 88 176 352 ...
                   13 26 52 104 208 15 30 60 120 240]);
    idx = find(Kb * Z_set >= K, 1);
    Z = Z_set(idx);
    K_ldpc = Kb * Z;
    parity_len = nb*Z - K_ldpc;
    expected_coded = K + parity_len;
    assert(length(coded) == expected_coded, ...
        sprintf('LDPC coded length: got %d, expected %d', length(coded), expected_coded));
    
    % Noiseless decode: create perfect LLRs
    llr_in = (1 - 2*coded) * 20;  % Strong LLRs
    [decoded, success] = ldpc_decode_block(llr_in, 8, K);
    assert(success, 'LDPC noiseless decode did not converge');
    assert(isequal(decoded(1:K), info), 'LDPC noiseless roundtrip failed');
    
    fprintf('  LDPC encode length: PASS\n');
    fprintf('  LDPC noiseless syndrome check: PASS\n');
    fprintf('  LDPC noiseless roundtrip: PASS\n');
end

%% ==================== OFDM Tests ====================
function test_ofdm()
    fprintf('--- OFDM Tests ---\n');
    
    NFFT = 64; NCP = 4;
    
    % Generate random QPSK symbols
    bits = randi([0 1], 80, 1);
    symbols = gray_map(bits, 2);  % QPSK
    
    [tx_sig, sc_idx, num_sym] = ofdm_modulate(symbols, NFFT, NCP);
    
    % Check signal length
    expected_len = num_sym * (NFFT + NCP);
    assert(length(tx_sig) == expected_len, 'OFDM tx signal length wrong');
    
    % Check subcarrier indices are valid
    assert(all(sc_idx >= 1 & sc_idx <= NFFT), 'Subcarrier indices out of range');
    
    % Noiseless channel test (flat fading h=[1])
    sigma2_w = 1e-10;
    h_taps = 1;
    rx_sig = tx_sig;  % No channel
    [eq_sym, snr_out] = ofdm_demod_equalize(rx_sig, h_taps, NFFT, NCP, num_sym, sc_idx, sigma2_w, 'zf');
    
    eq_sym = eq_sym(1:length(symbols));
    err = norm(eq_sym - symbols) / norm(symbols);
    assert(err < 1e-6, sprintf('OFDM noiseless roundtrip error: %.2e', err));
    
    fprintf('  OFDM modulate signal length: PASS\n');
    fprintf('  OFDM subcarrier indices: PASS\n');
    fprintf('  OFDM noiseless roundtrip (ZF): PASS\n');
end

%% ==================== Equalizer Tests ====================
function test_equalizers()
    fprintf('--- Equalizer Tests ---\n');
    
    NFFT = 64; NCP = 4;
    bits = randi([0 1], 60, 1);
    symbols = gray_map(bits, 2);
    [tx_sig, sc_idx, num_sym] = ofdm_modulate(symbols, NFFT, NCP);
    
    % Multipath channel with very low noise
    sigma2_w = 1e-8;
    L = 4;
    h = (randn(L,1) + 1j*randn(L,1))/sqrt(2*L);
    
    % Manual convolution
    hp = [h; zeros(length(tx_sig)-L, 1)];
    rx_sig = conv(tx_sig, h);
    rx_sig = rx_sig(1:length(tx_sig));
    
    % Test ZF
    [eq_zf, snr_zf] = ofdm_demod_equalize(rx_sig, h, NFFT, NCP, num_sym, sc_idx, sigma2_w, 'zf');
    eq_zf = eq_zf(1:length(symbols));
    
    % Test MMSE
    [eq_mmse, snr_mmse] = ofdm_demod_equalize(rx_sig, h, NFFT, NCP, num_sym, sc_idx, sigma2_w, 'mmse');
    eq_mmse = eq_mmse(1:length(symbols));
    
    % At very low noise, both should recover symbols well
    err_zf = norm(eq_zf - symbols) / norm(symbols);
    err_mmse = norm(eq_mmse - symbols) / norm(symbols);
    assert(err_zf < 0.1, sprintf('ZF eq error too high: %.4f', err_zf));
    assert(err_mmse < 0.1, sprintf('MMSE eq error too high: %.4f', err_mmse));
    
    % MMSE SNR should always be >= ZF SNR
    assert(all(snr_mmse(1:length(symbols)) >= snr_zf(1:length(symbols)) - 1e-6), ...
        'MMSE SNR should be >= ZF SNR');
    
    fprintf('  ZF equalizer low-noise recovery: PASS\n');
    fprintf('  MMSE equalizer low-noise recovery: PASS\n');
    fprintf('  MMSE SNR >= ZF SNR: PASS\n');
end

%% ==================== Rate Matching Tests ====================
function test_rate_matching()
    fprintf('--- Rate Matching Tests ---\n');
    
    coded = randi([0 1], 2060, 1);
    
    % No puncturing (E = coded_len)
    rm = rate_match(coded, 2060);
    assert(length(rm) == 2060, 'Rate match same length failed');
    assert(isequal(rm, coded), 'Rate match identity failed');
    
    % Repetition (E > coded_len)
    rm2 = rate_match(coded, 4000);
    assert(length(rm2) == 4000, 'Rate match repetition length failed');
    
    % Recovery
    llr_in = randn(2000, 1);
    llr_out = rate_recover(llr_in, 2060, 2000);
    assert(length(llr_out) == 2060, 'Rate recover length failed');
    assert(isequal(llr_out(1:2000), llr_in), 'Rate recover content mismatch');
    assert(all(llr_out(2001:end) == 0), 'Rate recover padding should be 0');
    
    fprintf('  Rate match identity: PASS\n');
    fprintf('  Rate match repetition: PASS\n');
    fprintf('  Rate recover: PASS\n');
end

%% ==================== Scrambling Tests ====================
function test_scrambling()
    fprintf('--- Scrambling Tests ---\n');
    
    len = 1000;
    seq = scramble_seq(len, 1, 0);
    assert(length(seq) == len, 'Scramble seq length wrong');
    assert(all(seq == 0 | seq == 1), 'Scramble seq not binary');
    
    % Same seed -> same sequence
    seq2 = scramble_seq(len, 1, 0);
    assert(isequal(seq, seq2), 'Scramble seq not deterministic');
    
    % Different RNTI -> different sequence
    seq3 = scramble_seq(len, 2, 0);
    assert(~isequal(seq, seq3), 'Different RNTI should give different seq');
    
    % Descramble LLR: scrambling then descrambling should preserve sign
    bits = randi([0 1], len, 1);
    scr_bits = xor(bits, seq);
    llr = 1 - 2*double(scr_bits);  % +1 for 0, -1 for 1
    llr_desc = descramble_llr(llr, seq);
    hard = double(llr_desc < 0);
    assert(isequal(hard, bits), 'Descramble LLR roundtrip failed');
    
    fprintf('  Scramble sequence length: PASS\n');
    fprintf('  Scramble deterministic: PASS\n');
    fprintf('  Scramble/Descramble roundtrip: PASS\n');
end
