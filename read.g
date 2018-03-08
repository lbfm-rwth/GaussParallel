# Loads all necessary functions into standard GAP, assuming access to 
# packages "IO" and "GAUSS".


LoadPackage("GAUSS");
LoadPackage("IO");
Read("./utils.g");
Read("./subprograms.g");
Read("./main_seq_trafo.g");
Read("./timing.g");
