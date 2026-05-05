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
