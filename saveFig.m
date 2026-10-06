function saveFig(h, file)
%SAVEFIG Export a single-column (3.4 in wide) figure as a vector PDF.
ht = get(h, 'UserData'); if isempty(ht) || ~isnumeric(ht), ht = 1.95; end
set(h, 'Units', 'inches', 'Position', [1 1 3.4 ht]);
if exist('exportgraphics', 'file')
    exportgraphics(h, file, 'ContentType', 'vector');
else
    set(h, 'PaperUnits', 'inches', 'PaperSize', [3.4 ht], 'PaperPosition', [0 0 3.4 ht]);
    print(h, file, '-dpdf');
end
end
