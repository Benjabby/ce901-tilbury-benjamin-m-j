function labelScatter(totals, metricX, metricY)
    cols = reshape([totals.col],3,[])';
    scatter([totals.(metricX)],[totals.(metricY)],46,cols,'filled');
    labelpoints([totals.(metricX)],[totals.(metricY)],[totals.name],'N',0.15,1);
end

