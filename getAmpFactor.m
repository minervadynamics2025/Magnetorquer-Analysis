function [ AmpFactor ] = getAmpFactor( r_mag, l_mag)

t_foot = linspace(0,0.15, 100)

AmpFactor =(l_mag * sqrt(r_mag^2 + (l_mag + t_foot).^2).*sqrt(r_mag^2 + t_foot.^2) );
AmpFactor =  AmpFactor./(((l_mag + t_foot) .* sqrt(r_mag^2 + t_foot.^2) - t_foot .* sqrt(r_mag^2 + (l_mag + t_foot).^2))* sqrt(r_mag^2 + 0.25 * l_mag^2))

figure;
plot(t_foot,AmpFactor)
end
