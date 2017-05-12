#include <faunus/faunus.h>
#include <faunus/ewald.h>
using namespace Faunus;
using namespace Faunus::Potential;

typedef Space<Geometry::Cuboid,PointParticle> Tspace;
typedef CombinedPairPotential<Coulomb, HardSphere> Tpairpot;

int main() {
  InputMap mcp("excess.json");                  // open user input file
  Tspace spc(mcp);                              // simulation space

  auto pot = Energy::Nonbonded<Tspace,Tpairpot>(mcp);

  pot.setSpace(spc);                            // share space w. hamiltonian

  spc.load("state",Tspace::RESIZE);             // load old config. from disk (if any)

  Analysis::CombinedAnalysis analyze(mcp,pot,spc);
  Move::Propagator<Tspace> mv(mcp,pot,spc);

  cout << atom.info() + spc.info() + pot.info() + "\n";

  MCLoop loop(mcp);                             // class for handling mc loops
  while ( loop[0] ) {
    while ( loop[1] ) {
      mv.move();                                // move!
      analyze.sample();
    }                                           // end of micro loop
    cout << loop.timing();
  }                                             // end of macro loop

  cout << loop.info() + mv.info() + analyze.info();
}
