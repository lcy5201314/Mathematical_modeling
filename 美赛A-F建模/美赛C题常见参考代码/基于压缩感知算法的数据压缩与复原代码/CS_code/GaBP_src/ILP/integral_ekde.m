function [out] = integral_ekde(in, weight,single)
    assert(weight ~= 0);
    if (nargin < 3)
        single = 0;
    end
    verify_kde(in,single);
    vij = -weight^2.*unique(getBW(in)).^2;
    if (abs(vij) < 1e-70)
%        assert(abs(vij) > 1e-70);
    end
    mij = -weight*getPoints(in);
    out = kde(mij./vij, 1./sqrt(vij));
    verify_kde(out,single);
end