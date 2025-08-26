% Backward compatibility wrapper
% The main driver script has been renamed to SchizoVis.m
% This thin wrapper simply invokes the new script so legacy calls to SchizoViz still work.
fprintf('[Compat] SchizoViz.m -> executing SchizoVis.m\n');
run('SchizoVis.m');