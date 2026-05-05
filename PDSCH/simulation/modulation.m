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
