#include "util/util.h"
#include "cudd/cudd.h"

namespace cudd {

class cudd_zdd 
{
public:
  cudd_zdd() 
  {
    ddmanager = Cudd_Init( 0, 0, CUDD_UNIQUE_SLOTS, CUDD_CACHE_SLOTS, 0 );
  }

  ~cudd_zdd()
  {
    Cudd_Quit( ddmanager );
  }

  void make_node()
  {
    DdNode *dd = Cudd_bddNewVar( ddmanager ); /*Create a new BDD variable*/
    Cudd_Ref( dd );

    //Cudd_PrintInfo( ddmanager, );

    int n = 1;
    int pr = 10;
    Cudd_PrintDebug( ddmanager, dd, n, pr );
  }

private:
  DdManager* ddmanager;
};

} // namespace cudd