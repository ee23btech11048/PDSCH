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
