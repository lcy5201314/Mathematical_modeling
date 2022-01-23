Usage:
DCVS(SeqName, NF, GOPSize, MKP, MCSP, FrameWidth, FrameHeight, Iteration)

SeqName: the name of test video sequence
NF: the number of test frames
GOPSize: GOP size
MKP: the measurement rate for each key frame
MCSP: the measurement rate for each CS frame
FrameWidth: frame width
FrameHeight: frame height
Iteration: the number of maximum iteration performed

Ex: 
DCVS('Coastguard_cif', 7, 3, 0.7, 0.4, 352, 288, 100);

The reconstructed frames will be saved as
RecCoastguard_cif000.bmp
RecCoastguard_cif001.bmp
.
.
.
RecCoastguard_cif006.bmp

