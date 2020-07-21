#include "cplusplus/cuddObj.hh"
#include "cudd/cuddInt.h"
#include <vector>
#include <string>
#include <iostream>
#include <unordered_map>
#include <array>
#include <utility>

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
    assert( Cudd_DebugCheck( cudd.getManager() ) == 0 );
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
    assert( Cudd_DebugCheck( cudd.getManager() ) == 0 );
    return node;
  }

  void deref_rec( ZDD const& node )
  {
    Cudd_RecursiveDerefZdd( cudd.getManager(), node.getNode() );
    assert( Cudd_DebugCheck( cudd.getManager() ) == 0 );
  }

  /* only decrease ref count of `node` but not recursively on its children */
  void deref( ZDD const& node )
  {
    deref( node.getNode() );
    assert( Cudd_DebugCheck( cudd.getManager() ) == 0 );
  }

  /* Create a node at level `var`, whose "then" child points to `T` and "else" child points to `E`,
   * or return the node if it already exists.
   * The result ZDD object is pass by copy, hence ref count is incresed by the copy constructor.
   * 
   * `T` and `E` should be lower than `var` (larger indices are lower)
   * Note the order of `T` and `E` is the same as in `bill`, but different from CUDD.
   */
  ZDD unique( uint32_t var, ZDD const& E, ZDD const& T )
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
    assert( ( f == empty || f == base || f.NodeReadIndex() < num_variables) && ( g == empty || g == base || g.NodeReadIndex() < num_variables ) );
    return ZDD( cudd, deref( join( f.getNode(), g.getNode() ) ) );
  }

  /* resulting sets are elements in f, but not superset of any element in g */
  /* \forall A \in result, A \in f and \forall B \in g, A \notsuperset B */
  ZDD nonsupersets( ZDD const& f, ZDD const& g )
  {
    assert( ( f == empty || f == base || f.NodeReadIndex() < num_variables) && ( g == empty || g == base || g.NodeReadIndex() < num_variables ) );
    return ZDD( cudd, deref( nonsupersets( f.getNode(), g.getNode() ) ) );
  }

  ZDD choose( ZDD const& f, uint32_t k )
  {
    assert( f == empty || f == base || f.NodeReadIndex() < num_variables );
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
  DdNode* deref( DdNode* node )
  {
    Cudd_Deref( node );
    return node;
    //cuddSatDec( node->ref );
    //if ( node->ref == 0 )
    //{
    //  DdManager * table = cudd.getManager();
    //  table->deadZ++;
    //  #ifdef DD_STATS
    //    table->nodesDropped++;
    //  #endif
    //  #ifdef DD_DEBUG
    //    assert(!cuddIsConstant(node));
    //  #endif
    //  int ord = table->permZ[node->index];
    //  table->subtableZ[ord].dead++;
    //  if (!cuddIsConstant(cuddT(node))) Cudd_RecursiveDerefZdd(table, cuddT(node));
    //  if (!cuddIsConstant(cuddE(node))) Cudd_RecursiveDerefZdd(table, cuddE(node));
    //}
    //return node; 
  }

  DdNode* unique( uint32_t var, DdNode* lo, DdNode* hi )
  { return ref( cuddZddGetNode( cudd.getManager(), var, hi, lo ) ); }

  DdNode* union_( DdNode* f, DdNode* g )
  { return ref( Cudd_zddUnion( cudd.getManager(), f, g ) ); }

  DdNode* intersection( DdNode* f, DdNode* g )
  { return ref( Cudd_zddIntersect( cudd.getManager(), f, g ) ); }

  DdNode* difference( DdNode* f, DdNode* g )
  { return ref( Cudd_zddDiff( cudd.getManager(), f, g ) ); }

#if 0
  DdNode* join_( DdManager * mgr, DdNode* f, DdNode* g ) { (void)mgr; return join( f, g ); }
  DdNode* nonsupersets_( DdManager * mgr, DdNode* f, DdNode* g ) { (void)mgr; return nonsupersets( f, g ); }
  DdNode* choose_( DdManager * mgr, DdNode* f, DdNode* g ) { (void)mgr; return choose( f, (uint64_t)g ); }
#else
  uint64_t join_ = 2;
  uint64_t nonsupersets_ = 4;
  uint64_t choose_ = 6;
#endif

  DdNode* join( DdNode* f, DdNode* g )
  {
    //std::cout<<"[i] compute join\n";
    /* terminal cases */
    DdNode* e = empty.getNode();
    DdNode* b = base.getNode();
    if ( f == e || g == e ) { return ref( e ); }
    if ( f == b ) { return ref( g ); }
    if ( g == b ) { return ref( f ); }

    /* cache lookup */
    //constexpr operations op = operations::zdd_join;
    //const auto it = computed_tables.at( op ).find( operand_t( f, g ) );
    //if ( it != computed_tables.at( op ).end() )
    //{
    //  if ( it->second->ref == 0 )
    //    { std::cout<<"    reclaim\n"; cuddReclaimZdd( cudd.getManager(), it->second ); }
    //  return ref( it->second );
    //}
    auto res = cache_lookup( cudd.getManager(), join_, f, g );
    if ( res != NULL ) { return ref( res ); }

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
    
    //computed_tables.at( op )[operand_t( f, g )] = r;
    cache_insert( cudd.getManager(), join_, f, g, r );
    assert( Cudd_DebugCheck( cudd.getManager() ) == 0 );
    return r;
  }

  DdNode* nonsupersets( DdNode* f, DdNode* g )
  {
    //std::cout<<"[i] compute nonsupersets\n";
    /* terminal cases */
    DdNode* e = empty.getNode();
    DdNode* b = base.getNode();
    if ( g == e ) { return ref( f ); }
    if ( f == e || g == b || f == g ) { return ref( e ); }

    if ( Cudd_NodeReadIndex( f ) > Cudd_NodeReadIndex( g ) )
      { return nonsupersets( f, cuddE( g ) ); }

    /* cache lookup */
    //constexpr operations op = operations::zdd_nonsupersets;
    //const auto it = computed_tables.at( op ).find( operand_t( f, g ) );
    //if ( it != computed_tables.at( op ).end() )
    //{
    //  if ( it->second->ref == 0 )
    //    { std::cout<<"    reclaim\n"; cuddReclaimZdd( cudd.getManager(), it->second ); }
    //  return ref( it->second );
    //}
    auto res = cache_lookup( cudd.getManager(), nonsupersets_, f, g );
    if ( res != NULL ) { return ref( res ); }

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
    
    //computed_tables.at( op )[operand_t( f, g )] = r;
    cache_insert( cudd.getManager(), nonsupersets_, f, g, r );
    assert( Cudd_DebugCheck( cudd.getManager() ) == 0 );
    return r;
  }

  DdNode* choose( DdNode* f, uint64_t k )
  {
    //std::cout<<"[i] compute choose\n";
    if ( Cudd_NodeReadIndex( f ) >= num_variables )
    {
      return ref( k > 0 ? empty.getNode() : base.getNode() );
    }
    if ( k == 1 ) { return ref( f ); }

    /* cache lookup */
    //constexpr operations op = operations::zdd_choose;
    //const auto it = computed_tables.at( op ).find( operand_t( f, (DdNode*)k ) );
    //if ( it != computed_tables.at( op ).end() )
    //{
    //  if ( it->second->ref == 0 )
    //    { std::cout<<"    reclaim\n"; cuddReclaimZdd( cudd.getManager(), it->second ); }
    //  return ref( it->second );
    //}
    auto res = cache_lookup( cudd.getManager(), choose_ + k * 2, f, f );
    if ( res != NULL ) { return ref( res ); }

    DdNode* n = choose( cuddE( f ), k ); /* don't take this var */
    if ( k > 0 )
    {
      DdNode* tmp = choose( cuddE( f ), k - 1 ); /* take this var */
      auto r = unique( Cudd_NodeReadIndex( f ), n, tmp );
      deref( n ); deref( tmp );
      //computed_tables.at( op )[operand_t( f, (DdNode*)k )] = r;
      cache_insert( cudd.getManager(), choose_ + k * 2, f, f, r );
      assert( Cudd_DebugCheck( cudd.getManager() ) == 0 );
      return r;
    }
    else /* k == 0 */
    {
      //computed_tables.at( op )[operand_t( f, (DdNode*)k )] = n;
      cache_insert( cudd.getManager(), choose_ + k * 2, f, f, n );
      assert( Cudd_DebugCheck( cudd.getManager() ) == 0 );
      return n;
    }
  }

private: /* operation cache */

#if 0 /* use our own cache */
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
  struct op_hash
  {
    std::size_t operator() ( operand_t const& ops ) const
    {
      return (std::size_t)ops.first + (std::size_t)ops.second;
    }
  };
  std::array<std::unordered_map<operand_t, DdNode*, op_hash>, operations::num_operations> computed_tables;

#else /* use CUDD's cache */
  DdNode * cache_lookup(
  DdManager * table,
  uint64_t op,
  DdNode * f,
  DdNode * g)
  {
    int posn;
    DdCache *en,*cache;
    DdNode *data;

    cache = table->cache;
  #ifdef DD_DEBUG
    if (cache == NULL) {
      return(NULL);
    }
  #endif

    posn = ddCHash2(op,f,g,table->cacheShift);
    en = &cache[posn];
    if (en->data != NULL && en->f==f && en->g==g && en->h==(ptruint)op) {
      data = Cudd_Regular(en->data);
      table->cacheHits++;
      if (data->ref == 0) {
        cuddReclaimZdd(table,data);
      }
      return(en->data);
    }

    /* Cache miss: decide whether to resize. */
    table->cacheMisses++;

    if (table->cacheSlack >= 0 && table->cacheHits > table->cacheMisses * table->minHit) {
      cuddCacheResize(table);
    }

    return(NULL);
  }

  void cache_insert(
  DdManager * table,
  uint64_t op,
  DdNode * f,
  DdNode * g,
  DdNode * data)
  {
    int posn;
    DdCache *entry;

    posn = ddCHash2(op,f,g,table->cacheShift);
    entry = &table->cache[posn];

    if (entry->data != NULL) {
      table->cachecollisions++;
    }
    table->cacheinserts++;

    entry->f = f;
    entry->g = g;
    entry->h = (ptruint) op;
    entry->data = data;
  #ifdef DD_CACHE_PROFILE
    entry->count++;
  #endif
  } 

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