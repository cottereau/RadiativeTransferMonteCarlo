function cells = lambMediumCells(I)
%LAMBMEDIUMCELLS  Body wave types of each medium of a lamb.Interface, with their speeds.
%   cells(m).types  {'P1' 'SV1' 'SH1'} | {'Pf'} | {'P2' 'SV2' 'SH2'};  cells(m).c  speeds;
%   cells(m).thc{b}  critical angles of type b (sin = c_b / every faster speed of both media).
md = {I.lower, I.upper};  nm = 1 + ~isempty(I.upper);  allc = [];
for m = 1:nm
    if md{m}.isFluid
        cells(m).types = {'Pf'};  cells(m).c = md{m}.alpha; %#ok<AGROW>
    else
        s = num2str(m);  cells(m).types = {['P' s], ['SV' s], ['SH' s]}; %#ok<AGROW>
        cells(m).c = [md{m}.alpha, md{m}.beta, md{m}.beta]; %#ok<AGROW>
    end
    allc = [allc, cells(m).c]; %#ok<AGROW>
end
for m = 1:nm
    for b = 1:numel(cells(m).types)
        cb = cells(m).c(b);  cells(m).thc{b} = asin(cb./allc(allc > cb*(1 + 1e-12)));
    end
end
end
