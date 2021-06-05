
systems = containers.Map;
systems('DCP') = @dehaze_he;
systems('Tarel') = @dehaze_tarel;

% blah blah
knowns.K = calib.P_rect{3};

% Infer the atmospheric light from transmission and final image
if predA==0
    predA = mean((predImage.*predT - img) ./ (predT-1),[1 2]);
end