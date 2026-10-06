function mu_c = constMuOrInf( B, H_0, N_d, Bmax )
%CONSTMUORINF Constant-mu_r fit, or Inf when the data exceed the mu_r -> inf limit.
if any(B(:) >= Bmax(:))
    mu_c = Inf;
else
    mu_c = constantPermeabilityFit(B, H_0, N_d);
end
end
