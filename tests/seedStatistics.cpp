//mst dependencies
#include "mstsystem.h"
#include "mstoptions.h"

//tpd dependencies
#include "benchmarkutilities.h"
#include "secondarystructure.h"
#include "utilities.h"

int main(int argc, char *argv[]) {
  MstOptions op;
  op.setTitle("Load seeds from a binary file, writes their statistics");
  op.addOption("pdb", "The peptide-protein complex that the seeds cover",true);
  op.addOption("p_id", "The ID of the peptide chain",true);
  op.addOption("bin_in", "Binary file with the current seeds",true);
  op.addOption("out", "The prefix for the new binary/csv files",true);
  op.setOptions(argc, argv);
  
  int seed_num = 100000;
  
  // Variables provided by user
  string bin_in = op.getString("bin_in");
  string out_prefix = op.getString("out");
  string pdb_path = op.getString("pdb");
  string p_id = op.getString("p_id");
  
  //// Make directories
  bool makeParents = true;
  // Fragment output folder
  string outDir = "termextension_output/";
  string covDir = "coverage_output/";
//  MstSys::cmkdir(outDir,makeParents);
//  MstSys::cmkdir(covDir,makeParents);
  
  Structure S(pdb_path);
  
  //report statistics
  seedStatistics stats(S, p_id, bin_in);
  stats.writeStatisticstoFile(outDir, out_prefix, seed_num);
  
  return 0;
}
