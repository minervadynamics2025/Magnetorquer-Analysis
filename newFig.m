function h = newFig( height )
%NEWFIG Single-column figure (3.4 in wide) with a print-friendly default style.
h = figure('Units', 'inches', 'Position', [1 1 3.4 height], 'Color', 'w');
set(h, 'UserData', height);
set(h, 'DefaultAxesFontName', 'Times New Roman', 'DefaultTextFontName', 'Times New Roman', ...
    'DefaultAxesFontSize', 7.5, 'DefaultTextFontSize', 7.5);
end
