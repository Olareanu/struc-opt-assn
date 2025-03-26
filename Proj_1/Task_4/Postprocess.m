function [] = Postprocess(fem,opts)
    viz2D3D_line_deformed(fem,opts,1,'Mag'); 
    viz2D3D_line_stresses(fem,opts,1);
end