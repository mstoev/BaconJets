LIBRARY := SUHH2BaconJets
#DICT := include/ReconstructionHypothesis.h include/SUHH2common_LinkDef.h
#TEST := 1
USERLDFLAGS := -lSUHH2core -lSUHH2bacondataformats -lSUHH2common -lSUHH2JetMETObjects
PAR := 1
include ../Makefile.common

