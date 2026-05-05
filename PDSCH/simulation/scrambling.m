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
