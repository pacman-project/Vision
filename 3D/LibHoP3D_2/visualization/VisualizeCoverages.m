function [] = VisualizeCoverages(CoverageOut1, CoverageOut2, CoverageOut3, coverageOverall, len) 

    CoverageOut1 = CoverageOut1(1:len);
    CoverageOut2 = CoverageOut2(1:len);
    CoverageOut3 = CoverageOut3(1:len);
    
    CoverageOut1 = CoverageOut1 / coverageOverall;
    CoverageOut2 = CoverageOut2 / coverageOverall;
    CoverageOut3 = CoverageOut3 / coverageOverall;
    
    CoverageOut1 = cumsum(CoverageOut1);
    CoverageOut2 = cumsum(CoverageOut2);
    CoverageOut3 = cumsum(CoverageOut3);
    
    x = 1:len;
    
    plot(x, CoverageOut1, 'blue');
    hold on
    plot(x, CoverageOut2, 'green');
    plot(x, CoverageOut3, 'red');
    hold off

end