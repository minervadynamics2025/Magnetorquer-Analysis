function styleMarkers( h, c )
%STYLEMARKERS Open circular markers with thin, cap-less error bars in colour c.
set(h, 'MarkerSize', 4, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', c, 'Color', c, 'LineWidth', 1.0);
try, set(h, 'CapSize', 0); catch, end     % CapSize exists in MATLAB R2016b+
end
