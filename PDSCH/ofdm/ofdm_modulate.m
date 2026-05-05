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
