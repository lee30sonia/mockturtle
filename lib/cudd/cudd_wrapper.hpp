#include "cplusplus/cuddObj.hh"
#include "cudd/cuddInt.h"
#include <vector>
#include <string>
#include <iostream>
#include <unordered_map>
#include <array>

namespace cudd {

class cudd_zdd 
{
public:
  cudd_zdd( uint32_t num_variables ) 
    : num_variables( num_variables ), empty( cudd.zddZero() ), base( cudd.zddOne( INT_MAX ) )
  {
    elementaries.reserve( num_variables );
    for ( auto i = 0u; i < num_variables; ++i )
    {
      elementaries.emplace_back( base.Change( i ) );
    }
    tautologies.reserve( num_variables );
    for ( auto i = 0u; i < num_variables; ++i )
    {
      tautologies.emplace_back( cudd.zddOne( i ) );
    }
    assert( cudd.ReadZddSize() == int( num_variables ) );
    //cudd.makeVerbose();
  }

  ~cudd_zdd()
  {
    cudd.makeTerse(); 
  }

  void print( ZDD const& node, std::string const& name = "", int verbosity = 4 )
  {
    std::cout << name << " (" << node.getNode() << ", level = " << Cudd_NodeReadIndex( node.getNode() ) << ")";
    node.print( num_variables, verbosity );
  }

  ZDD ref( ZDD const& node )
  {
    Cudd_Ref( node.getNode() );
    return node;
  }

  void deref_rec( ZDD const& node )
  {
    Cudd_RecursiveDerefZdd( cudd.getManager(), node.getNode() );
  }

  /* only decrease ref count of `node` but not recursively on its children */
  void deref( ZDD const& node )
  {
    Cudd_Deref( node.getNode() );
  }

  /* Create a node at level `var`, whose "then" child points to `T` and "else" child points to `E`,
   * or return the node if it already exists.
   * The result ZDD object is pass by copy, hence ref count is incresed by the copy constructor.s
   * 
   * `T` and `E` should be lower than `var` (larger indices are lower)
   */
  ZDD create( uint32_t var, ZDD const& T, ZDD const& E )
  {
    assert( T.NodeReadIndex() > var && E.NodeReadIndex() > var );
    return ZDD( cudd, cuddZddGetNode( cudd.getManager(), var, T.getNode(), E.getNode() ) );
  }

public: /* basic getters; does NOT increase ref count */
  ZDD& bottom() { return empty; }
  ZDD& top() { return base; }

  ZDD& elementary( uint32_t var )
  {
    assert( var < num_variables );
    return elementaries[var];
  }

  ZDD& tautology( uint32_t var = 0 )
  {
    assert( var < num_variables );
    return tautologies[var];
  }

public: /* operations provided by CUDD, wrapped with `bill` function names */
  ZDD union_( ZDD const& f, ZDD const& g )
  {
    return f.Union( g );
  }

  ZDD intersection( ZDD const& f, ZDD const& g )
  {
    return f.Intersect( g );
  }

  ZDD difference( ZDD const& f, ZDD const& g )
  {
    return f.Diff( g );
  }

public: /* operations provided by CUDD */
  /* union every pair of subsets in f and g */
  ZDD join( ZDD const& f, ZDD const& g )
  {
    assert( f.NodeReadIndex() < num_variables && g.NodeReadIndex() < num_variables );
    return ZDD( cudd, deref( join( f.getNode(), g.getNode() ) ) );
  }

  /* resulting sets are elements in f, but not superset of any element in g */
  /* \forall A \in result, A \in f and \forall B \in g, A \notsuperset B */
  ZDD nonsupersets( ZDD const& f, ZDD const& g )
  {
    assert( f.NodeReadIndex() < num_variables && g.NodeReadIndex() < num_variables );
    return ZDD( cudd, deref( nonsupersets( f.getNode(), g.getNode() ) ) );
  }

  ZDD choose( ZDD const& f, uint32_t k )
  {
    assert( f.NodeReadIndex() < num_variables );
    return ZDD( cudd, deref( choose( f.getNode(), k ) ) );
  }

  //edivide
  //maximal
  //meet
  //nonsubsets

private: /* implementation details */
  DdNode* lo( DdNode* f ) { return cuddE( f ); }
  DdNode* hi( DdNode* f ) { return cuddT( f ); }
  DdNode* ref( DdNode* node ) { Cudd_Ref( node ); return node; }
  DdNode* deref( DdNode* node ) { Cudd_Deref( node ); return node; }

  /* note the order is the same as in `bill` but different than `create` */
  DdNode* unique( uint32_t var, DdNode* lo, DdNode* hi )
  { return ref( cuddZddGetNode( cudd.getManager(), var, hi, lo ) ); }

  DdNode* union_( DdNode* f, DdNode* g )
  { return ref( Cudd_zddUnion( cudd.getManager(), f, g ) ); }

  DdNode* intersection( DdNode* f, DdNode* g )
  { return ref( Cudd_zddIntersect( cudd.getManager(), f, g ) ); }

  DdNode* difference( DdNode* f, DdNode* g )
  { return ref( Cudd_zddDiff( cudd.getManager(), f, g ) ); }

  DdNode* join( DdNode* f, DdNode* g )
  {
    /* terminal cases */
    DdNode* e = empty.getNode();
    DdNode* b = base.getNode();
    if ( f == e || g == e ) { return ref( e ); }
    if ( f == b ) { return ref( g ); }
    if ( g == b ) { return ref( f ); }

    uint32_t var = std::min( Cudd_NodeReadIndex( f ), Cudd_NodeReadIndex( g ) );
    DdNode* lo = 0;
    DdNode* hi = 0;
    if ( Cudd_NodeReadIndex( f ) < Cudd_NodeReadIndex( g ) )
    {
      lo = join( cuddE( f ), g );
      hi = join( cuddT( f ), g );
    }
    else if ( Cudd_NodeReadIndex( f ) > Cudd_NodeReadIndex( g ) )
    {
      lo = join( cuddE( g ), f );
      hi = join( cuddT( g ), f );
    }
    else /* f_var == g_var */
    {
      DdNode* tmp0 = union_( cuddE( g ), cuddT( g ) );
      DdNode* tmp1 = join( cuddT( f ), tmp0 );
      deref( tmp0 );
      DdNode* tmp2 = join( cuddE( f ), cuddT( g ) );
      hi = union_( tmp1, tmp2 );
      deref( tmp1 ); deref( tmp2 );
      lo = join( cuddE( f ), cuddE( g ) );
    }
    auto r = unique( var, lo, hi );
    deref( lo ); deref( hi );
    return r;
  }

  DdNode* nonsupersets( DdNode* f, DdNode* g )
  {
    /* terminal cases */
    DdNode* e = empty.getNode();
    DdNode* b = base.getNode();
    if ( g == e ) { return ref( f ); }
    if ( f == e || g == b || f == g ) { return ref( e ); }

    if ( Cudd_NodeReadIndex( f ) > Cudd_NodeReadIndex( g ) )
      { return nonsupersets( f, cuddE( g ) ); }

    /* cache lookup */
    /*constexpr operations op = operations::zdd_nonsupersets;
    const auto it = computed_tables.at( op ).find( hash( f, g ) );
    if ( it != computed_tables.at( op ).end() )
    {
      if ( it->second->ref == 0 )
        { cuddReclaimZdd( cudd.getManager(), it->second ); }
      return ref( it->second );
    }*/

    /* recursive computation */
    DdNode* lo = 0;
    DdNode* hi = 0;
    if ( Cudd_NodeReadIndex( f ) < Cudd_NodeReadIndex( g ) )
    {
      lo = nonsupersets( cuddE( f ), g );
      hi = nonsupersets( cuddT( f ), g );
    }
    else /* f_var == g_var */
    {
      DdNode* tmp0 = nonsupersets( cuddT( f ), cuddT( g ) );
      DdNode* tmp1 = nonsupersets( cuddT( f ), cuddE( g ) );
      hi = intersection( tmp0, tmp1 );
      deref( tmp0 ); deref( tmp1 );
      lo = nonsupersets( cuddE( f ), cuddE( g ) );
    }
    auto r = unique( Cudd_NodeReadIndex( f ), lo, hi );
    deref( lo ); deref( hi );
    //computed_tables.at( op )[hash( f, g )] = r;
    return r;
  }

  DdNode* choose( DdNode* f, uint32_t k )
  {
    if ( Cudd_NodeReadIndex( f ) >= num_variables )
    {
      return ref( k > 0 ? empty.getNode() : base.getNode() );
    }
    if ( k == 1 ) { return ref( f ); }

    DdNode* n = choose( cuddE( f ), k ); /* don't take this var */
    if ( k > 0 )
    {
      DdNode* tmp = choose( cuddE( f ), k - 1 ); /* take this var */
      auto r = unique( Cudd_NodeReadIndex( f ), n, tmp );
      deref( n ); deref( tmp );
      return r;
    }
    else /* k == 0 */
    {
      return n;
    }
  }

#if 0
private: /* operation cache */
  enum operations : uint32_t 
  {
    zdd_choose,
    //zdd_difference,
    //zdd_edivide,
    //zdd_intersection,
    zdd_join,
    //zdd_maximal,
    //zdd_meet,
    //zdd_nonsubsets,
    zdd_nonsupersets,
    //zdd_union,
    num_operations
  };
  using operand_t = std::pair<DdNode*, DdNode*>;
  std::array<std::unordered_map<operand_t, DdNode*>, operations::num_operations> computed_tables;
#endif

private:
  Cudd cudd; /* the CUDD manager */
  uint32_t num_variables;
  ZDD empty; /* the empty family {} (minterms: none; constant 0) */
  ZDD base; /* the unit family {{}} (minterms: the all-zero cube) */
  std::vector<ZDD> elementaries; /* the single-set family of the single-element set {{i}} (minterms: 0...010...0) */
  std::vector<ZDD> tautologies; /* every combinations of variables >= var (minterms: 0...0-...-) */
};

} // namespace cudd