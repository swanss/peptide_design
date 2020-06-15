#include "params.h"

RotamerLibrary* getRotamerLibrary(int argc, char *argv[]) {
	MstOptions opts;
	opts.addOption("rotLib", "Location of rotmaer library.", false);
	opts.setOptions(argc, argv);
	string rotLibFile = opts.getString("rotLib", "/home/grigoryanlab/library/MST/testfiles/rotlib.bin");
	return new RotamerLibrary(rotLibFile);
}
