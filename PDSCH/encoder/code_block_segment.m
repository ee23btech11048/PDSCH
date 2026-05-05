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
