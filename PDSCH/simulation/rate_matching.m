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
