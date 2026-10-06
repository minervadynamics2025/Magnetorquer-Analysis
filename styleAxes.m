function styleAxes( ax )
%STYLEAXES Recessive grid and axes: thin light grid, inward ticks, no top/right box.
INK2 = [82 81 78]/255;
set(ax, 'Box', 'off', 'TickDir', 'in', 'LineWidth', 0.6, 'XColor', INK2, 'YColor', INK2, ...
    'XGrid', 'on', 'YGrid', 'on', 'GridColor', [228 227 223]/255, 'GridAlpha', 1, 'Layer', 'bottom');
lg = findobj(ancestor(ax, 'figure'), 'Type', 'Legend');
set(lg, 'Box', 'off', 'FontSize', 7);
end
