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