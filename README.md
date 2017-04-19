# FSCW_winding_layout
Matlab M-script to analyze winding layout of fractional-slot concentrated winding permanent magnet machines.

Created in Matlab 2016b.
'cros.m' is a function where you input the stator slot number and rotor pole number, it will output the optimal winding layout for that machine according to Cros' method and display illustration of phasor alignment and winding layout. It will also tell you if it is impossible to build a reasonable machine with that combination of pole/slots.
'cros_batch.m' calls cros.m to evaluate 1200 different designs (execution: a few seconds) and display the key performance index 'winding factor' for all designs. It also gives illustrations of two reasonably designed machines.

Files are extensive commented, look inside for more hints about how they work.
