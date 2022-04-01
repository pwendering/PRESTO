function gurobiSense = translateSenseCobra2Gurobi(cobraSense)
    gurobiSense = strrep(cobraSense','E','=');
    gurobiSense = strrep(gurobiSense,'L','<');
    gurobiSense = strrep(gurobiSense,'G','>')';
end