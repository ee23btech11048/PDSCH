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
