
systems = containers.Map;
systems('DCP') = @dehaze_he;
systems('Tarel') = @dehaze_tarel;

% blah blah
knowns.K = calib.P_rect{3};

% use timeit for full function.
% timeTotal =/= timeImage + timeA