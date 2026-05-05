function decoded = viterbi_decode(llr, K, gen_poly)
    persistent cached_K cached_gen_poly cached_nxt cached_outp cached_ns cached_ng cached_mem
    
    llr = llr(:).';
    L = floor(length(llr) / length(gen_poly));
    
    if isempty(cached_K) || cached_K ~= K || ~isequal(cached_gen_poly, gen_poly)
        cached_K = K; 
        cached_gen_poly = gen_poly;
        cached_ng = length(gen_poly);
        cached_mem = K - 1;
        cached_ns = 2^cached_mem;
        
        gb = zeros(cached_ng, K);
        for g = 1:cached_ng
            % gen_poly values like 133, 171 are OCTAL representations stored as decimal
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
        
        cached_nxt = zeros(cached_ns, 2);
        cached_outp = zeros(cached_ns, 2, cached_ng);
        for s = 0:cached_ns-1
            sb = bitget(s, cached_mem:-1:1);
            for inp = 0:1
                reg = [inp sb];
                nsb = reg(1:cached_mem);
                cached_nxt(s+1, inp+1) = 0;
                for b = 1:cached_mem
                    cached_nxt(s+1, inp+1) = cached_nxt(s+1, inp+1) + nsb(b) * 2^(cached_mem-b);
                end
                for g = 1:cached_ng
                    cached_outp(s+1, inp+1, g) = mod(sum(reg .* gb(g,:)), 2);
                end
            end
        end
    end
    
    nxt = cached_nxt;
    outp = cached_outp;
    ns = cached_ns;
    ng = cached_ng;
    mem = cached_mem;

    pm = -inf(ns, 1);
    pm(1) = 0;
    surv = zeros(ns, L);
    sh = zeros(ns, L);
    
    for t = 1:L
        npm = -inf(ns, 1);
        nsu = zeros(ns, 1);
        nsh = zeros(ns, 1);
        bl = llr((t-1)*ng+1 : t*ng);
        for s = 0:ns-1
            if pm(s+1) == -inf, continue; end
            for inp = 0:1
                nst = nxt(s+1, inp+1);
                bm = 0;
                for g = 1:ng
                    ex = 1 - 2*outp(s+1, inp+1, g);
                    bm = bm + bl(g) * ex;
                end
                cand = pm(s+1) + bm;
                if cand > npm(nst+1)
                    npm(nst+1) = cand;
                    nsu(nst+1) = inp;
                    nsh(nst+1) = s;
                end
            end
        end
        pm = npm;
        surv(:,t) = nsu;
        sh(:,t) = nsh;
    end
    decoded = zeros(L, 1);
    s = 0;
    for t = L:-1:1
        decoded(t) = surv(s+1, t);
        s = sh(s+1, t);
    end
    decoded = decoded(1:L-mem);
end
