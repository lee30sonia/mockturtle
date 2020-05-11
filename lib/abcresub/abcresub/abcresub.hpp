/*!
  \file abcresub.hpp
  \brief Extracted from ABC
         https://github.com/berkeley-abc/abc

  \author Alan Mishchenko (UC Berkeley)
*/

#pragma once

#include "abc_global.hpp"
#include "vec_ptr.hpp"
#include "vec_int.hpp"
#include "vec_wrd.hpp"
#include "tt.hpp"

namespace abcresub
{

typedef struct Gia_ResbMan_t_ Gia_ResbMan_t;

struct Gia_ResbMan_t_
{
    int         nWords;
    int         nLimit;
    int         fDebug;
    int         fVerbose;
    Vec_Ptr_t * vDivs;
    Vec_Int_t * vGates;
    Vec_Int_t * vUnateLits[2];
    Vec_Int_t * vNotUnateVars[2];
    Vec_Int_t * vUnatePairs[2];
    Vec_Int_t * vBinateVars;
    Vec_Int_t * vUnateLitsW[2];
    Vec_Int_t * vUnatePairsW[2];
    word *      pSets[2];
    word *      pDivA;
    word *      pDivB;
    Vec_Wrd_t * vSims;
};

inline Gia_ResbMan_t * Gia_ResbAlloc( int nWords )
{
  Gia_ResbMan_t * p   = ABC_CALLOC( Gia_ResbMan_t, 1 );
  p->nWords           = nWords;
  p->vUnateLits[0]    = Vec_IntAlloc( 100 );
  p->vUnateLits[1]    = Vec_IntAlloc( 100 );
  p->vNotUnateVars[0] = Vec_IntAlloc( 100 );
  p->vNotUnateVars[1] = Vec_IntAlloc( 100 );
  p->vUnatePairs[0]   = Vec_IntAlloc( 100 );
  p->vUnatePairs[1]   = Vec_IntAlloc( 100 );
  p->vUnateLitsW[0]   = Vec_IntAlloc( 100 );
  p->vUnateLitsW[1]   = Vec_IntAlloc( 100 );
  p->vUnatePairsW[0]  = Vec_IntAlloc( 100 );
  p->vUnatePairsW[1]  = Vec_IntAlloc( 100 );
  p->vBinateVars      = Vec_IntAlloc( 100 );
  p->vGates           = Vec_IntAlloc( 100 );
  p->vDivs            = Vec_PtrAlloc( 100 );
  p->pSets[0]         = ABC_CALLOC( word, nWords );
  p->pSets[1]         = ABC_CALLOC( word, nWords );
  p->pDivA            = ABC_CALLOC( word, nWords );
  p->pDivB            = ABC_CALLOC( word, nWords );
  p->vSims            = Vec_WrdAlloc( 100 );
  return p;
}

inline void Gia_ResbInit( Gia_ResbMan_t * p, Vec_Ptr_t * vDivs, int nWords, int nLimit, int fDebug, int fVerbose )
{
  assert( p->nWords == nWords );
  p->nLimit   = nLimit;
  p->fDebug   = fDebug;
  p->fVerbose = fVerbose;
  Abc_TtCopy( p->pSets[0], (word *)Vec_PtrEntry(vDivs, 0), nWords, 0 );
  Abc_TtCopy( p->pSets[1], (word *)Vec_PtrEntry(vDivs, 1), nWords, 0 );
  Vec_PtrClear( p->vDivs );
  Vec_PtrAppend( p->vDivs, vDivs );
  Vec_IntClear( p->vGates );
  Vec_IntClear( p->vUnateLits[0]    );
  Vec_IntClear( p->vUnateLits[1]    );
  Vec_IntClear( p->vNotUnateVars[0] );
  Vec_IntClear( p->vNotUnateVars[1] );
  Vec_IntClear( p->vUnatePairs[0]   );
  Vec_IntClear( p->vUnatePairs[1]   );
  Vec_IntClear( p->vUnateLitsW[0]   );
  Vec_IntClear( p->vUnateLitsW[1]   );
  Vec_IntClear( p->vUnatePairsW[0]  );
  Vec_IntClear( p->vUnatePairsW[1]  );
  Vec_IntClear( p->vBinateVars      );
}

inline void Gia_ResbFree( Gia_ResbMan_t * p )
{
  Vec_IntFree( p->vUnateLits[0]    );
  Vec_IntFree( p->vUnateLits[1]    );
  Vec_IntFree( p->vNotUnateVars[0] );
  Vec_IntFree( p->vNotUnateVars[1] );
  Vec_IntFree( p->vUnatePairs[0]   );
  Vec_IntFree( p->vUnatePairs[1]   );
  Vec_IntFree( p->vUnateLitsW[0]   );
  Vec_IntFree( p->vUnateLitsW[1]   );
  Vec_IntFree( p->vUnatePairsW[0]  );
  Vec_IntFree( p->vUnatePairsW[1]  );
  Vec_IntFree( p->vBinateVars      );
  Vec_IntFree( p->vGates           );
  Vec_WrdFree( p->vSims            );
  Vec_PtrFree( p->vDivs            );
  ABC_FREE( p->pSets[0] );
  ABC_FREE( p->pSets[1] );
  ABC_FREE( p->pDivA );
  ABC_FREE( p->pDivB );
  ABC_FREE( p );
}

/**Function*************************************************************

  Synopsis    [Print resubstitution.]

  Description []
               
  SideEffects []

  SeeAlso     []

***********************************************************************/
inline void Gia_ManResubPrintNode( Vec_Int_t * vRes, int nVars, int Node, int fCompl )
{
    extern void Gia_ManResubPrintLit( Vec_Int_t * vRes, int nVars, int iLit );

    int iLit0 = Vec_IntEntry( vRes, 2*Node + 0 );
    int iLit1 = Vec_IntEntry( vRes, 2*Node + 1 );
    assert( iLit0 != iLit1 );
    if ( iLit0 > iLit1 && Abc_LitIsCompl(fCompl) ) // xor
    {
        printf( "~" );
        fCompl = 0;
    }
    printf( "(" );
    Gia_ManResubPrintLit( vRes, nVars, Abc_LitNotCond(iLit0, fCompl) );
    printf( " %c ", iLit0 > iLit1 ? '^' : (fCompl ? '|' : '&') );
    Gia_ManResubPrintLit( vRes, nVars, Abc_LitNotCond(iLit1, fCompl) );
    printf( ")" );
}

inline void Gia_ManResubPrintLit( Vec_Int_t * vRes, int nVars, int iLit )
{
    if ( Abc_Lit2Var(iLit) < nVars )
    {
        if ( nVars < 26 )
            printf( "%s%c", Abc_LitIsCompl(iLit) ? "~":"", 'a' + Abc_Lit2Var(iLit)-2 );
        else
            printf( "%si%d", Abc_LitIsCompl(iLit) ? "~":"", Abc_Lit2Var(iLit)-2 );
    }
    else
        Gia_ManResubPrintNode( vRes, nVars, Abc_Lit2Var(iLit) - nVars, Abc_LitIsCompl(iLit) );
}

inline int Gia_ManResubPrint( Vec_Int_t * vRes, int nVars )
{
    int iTopLit;
    if ( Vec_IntSize(vRes) == 0 )
        return printf( "none" );
    assert( Vec_IntSize(vRes) % 2 == 1 );
    iTopLit = Vec_IntEntryLast(vRes);
    if ( iTopLit == 0 )
        return printf( "const0" );
    if ( iTopLit == 1 )
        return printf( "const1" );
    Gia_ManResubPrintLit( vRes, nVars, iTopLit );
    return 0;
}

/**Function*************************************************************

  Synopsis    [Verify resubstitution.]

  Description []
               
  SideEffects []

  SeeAlso     []

***********************************************************************/
int Gia_ManResubVerify( Gia_ResbMan_t * p )
{
    int nVars = Vec_PtrSize(p->vDivs);
    int iTopLit, RetValue;
    word * pDivRes; 
    if ( Vec_IntSize(p->vGates) == 0 )
        return -1;
    iTopLit = Vec_IntEntryLast(p->vGates);
    assert( iTopLit >= 0 );
    if ( iTopLit == 0 )
        return Abc_TtIsConst0( p->pSets[1], p->nWords );
    if ( iTopLit == 1 )
        return Abc_TtIsConst0( p->pSets[0], p->nWords );
    if ( Abc_Lit2Var(iTopLit) < nVars )
    {
        assert( Vec_IntSize(p->vGates) == 1 );
        pDivRes = (word *)Vec_PtrEntry( p->vDivs, Abc_Lit2Var(iTopLit) );
    }
    else
    {
        int i, iLit0, iLit1;
        assert( Vec_IntSize(p->vGates) > 1 );
        assert( Vec_IntSize(p->vGates) % 2 == 1 );
        assert( Abc_Lit2Var(iTopLit)-nVars == Vec_IntSize(p->vGates)/2-1 );
        Vec_WrdFill( p->vSims, p->nWords * Vec_IntSize(p->vGates)/2, 0 );
        Vec_IntForEachEntryDouble( p->vGates, iLit0, iLit1, i )
        {
            int iVar0 = Abc_Lit2Var(iLit0);
            int iVar1 = Abc_Lit2Var(iLit1);
            word * pDiv0 = iVar0 < nVars ? (word *)Vec_PtrEntry(p->vDivs, iVar0) : Vec_WrdEntryP(p->vSims, p->nWords*(iVar0 - nVars));
            word * pDiv1 = iVar1 < nVars ? (word *)Vec_PtrEntry(p->vDivs, iVar1) : Vec_WrdEntryP(p->vSims, p->nWords*(iVar1 - nVars));
            word * pDiv  = Vec_WrdEntryP(p->vSims, p->nWords*i/2);
            if ( iVar0 < iVar1 )
                Abc_TtAndCompl( pDiv, pDiv0, Abc_LitIsCompl(iLit0), pDiv1, Abc_LitIsCompl(iLit1), p->nWords );
            else if ( iVar0 > iVar1 )
            {
                assert( !Abc_LitIsCompl(iLit0) );
                assert( !Abc_LitIsCompl(iLit1) );
                Abc_TtXor( pDiv, pDiv0, pDiv1, p->nWords, 0 );
            }
            else assert( 0 );
        }
        pDivRes = Vec_WrdEntryP( p->vSims, p->nWords*(Vec_IntSize(p->vGates)/2-1) );
    }
    if ( Abc_LitIsCompl(iTopLit) )
        RetValue = !Abc_TtIntersectOne(p->pSets[1], 0, pDivRes, 0, p->nWords) && !Abc_TtIntersectOne(p->pSets[0], 0, pDivRes, 1, p->nWords);
    else
        RetValue = !Abc_TtIntersectOne(p->pSets[0], 0, pDivRes, 0, p->nWords) && !Abc_TtIntersectOne(p->pSets[1], 0, pDivRes, 1, p->nWords);
    return RetValue;
}

#if 0
/**Function*************************************************************

  Synopsis    [Construct AIG manager from gates.]

  Description []
               
  SideEffects []

  SeeAlso     []

***********************************************************************/
inline int Gia_ManGetVar( Gia_Man_t * pNew, Vec_Int_t * vUsed, int iVar )
{
    if ( Vec_IntEntry(vUsed, iVar) == -1 )
        Vec_IntWriteEntry( vUsed, iVar, Gia_ManAppendCi(pNew) );
    return Vec_IntEntry(vUsed, iVar);
}

inline Gia_Man_t * Gia_ManConstructFromGates( Vec_Int_t * vGates, int nVars )
{
    int i, iLit0, iLit1, iTopLit, iLitRes;
    Gia_Man_t * pNew;
    if ( Vec_IntSize(vGates) == 0 )
        return NULL;
    assert( Vec_IntSize(vGates) % 2 == 1 );
    pNew = Gia_ManStart( 100 );
    pNew->pName = Abc_UtilStrsav( "resub" );
    iTopLit = Vec_IntEntryLast( vGates );
    if ( iTopLit == 0 || iTopLit == 1 )
        iLitRes = 0;
    else if ( Abc_Lit2Var(iTopLit) < nVars )
    {
        assert( Vec_IntSize(vGates) == 1 );
        iLitRes = Gia_ManAppendCi(pNew);
    }
    else
    {
        Vec_Int_t * vUsed = Vec_IntStartFull( nVars );
        Vec_Int_t * vCopy = Vec_IntAlloc( Vec_IntSize(vGates)/2 );
        assert( Vec_IntSize(vGates) > 1 );
        assert( Vec_IntSize(vGates) % 2 == 1 );
        assert( Abc_Lit2Var(iTopLit)-nVars == Vec_IntSize(vGates)/2-1 );
        Vec_IntForEachEntryDouble( vGates, iLit0, iLit1, i )
        {
            int iVar0 = Abc_Lit2Var(iLit0);
            int iVar1 = Abc_Lit2Var(iLit1);
            int iRes0 = iVar0 < nVars ? Gia_ManGetVar(pNew, vUsed, iVar0) : Vec_IntEntry(vCopy, iVar0 - nVars);
            int iRes1 = iVar1 < nVars ? Gia_ManGetVar(pNew, vUsed, iVar1) : Vec_IntEntry(vCopy, iVar1 - nVars);
            if ( iVar0 < iVar1 )
                iLitRes = Gia_ManAppendAnd( pNew, Abc_LitNotCond(iRes0, Abc_LitIsCompl(iLit0)), Abc_LitNotCond(iRes1, Abc_LitIsCompl(iLit1)) );
            else if ( iVar0 > iVar1 )
            {
                assert( !Abc_LitIsCompl(iLit0) );
                assert( !Abc_LitIsCompl(iLit1) );
                iLitRes = Gia_ManAppendXor( pNew, Abc_LitNotCond(iRes0, Abc_LitIsCompl(iLit0)), Abc_LitNotCond(iRes1, Abc_LitIsCompl(iLit1)) );
            }
            else assert( 0 );
            Vec_IntPush( vCopy, iLitRes );
        }
        assert( Vec_IntSize(vCopy) == Vec_IntSize(vGates)/2 );
        iLitRes = Vec_IntEntry( vCopy, Vec_IntSize(vGates)/2-1 );
        Vec_IntFree( vUsed );
        Vec_IntFree( vCopy );
    }
    Gia_ManAppendCo( pNew, Abc_LitNotCond(iLitRes, Abc_LitIsCompl(iTopLit)) );
    return pNew;
}
#endif

/**Function*************************************************************

  Synopsis    [Perform resubstitution.]

  Description []
               
  SideEffects []

  SeeAlso     []

***********************************************************************/
inline int Gia_ManFindFirstCommonLit( Vec_Int_t * vArr1, Vec_Int_t * vArr2, int fVerbose )
{
    int * pBeg1 = vArr1->pArray;
    int * pBeg2 = vArr2->pArray;
    int * pEnd1 = vArr1->pArray + vArr1->nSize;
    int * pEnd2 = vArr2->pArray + vArr2->nSize;
    int * pStart1 = vArr1->pArray;
    int * pStart2 = vArr2->pArray;
    int nRemoved = 0;
    while ( pBeg1 < pEnd1 && pBeg2 < pEnd2 )
    {
        if ( Abc_Lit2Var(*pBeg1) == Abc_Lit2Var(*pBeg2) )
        { 
            if ( *pBeg1 != *pBeg2 ) 
                return *pBeg1; 
            else
                pBeg1++, pBeg2++;
            nRemoved++;
        }
        else if ( *pBeg1 < *pBeg2 )
            *pStart1++ = *pBeg1++;
        else 
            *pStart2++ = *pBeg2++;
    }
    while ( pBeg1 < pEnd1 )
        *pStart1++ = *pBeg1++;
    while ( pBeg2 < pEnd2 )
        *pStart2++ = *pBeg2++;
    Vec_IntShrink( vArr1, pStart1 - vArr1->pArray );
    Vec_IntShrink( vArr2, pStart2 - vArr2->pArray );
    if ( fVerbose ) printf( "Removed %d duplicated entries.  Array1 = %d.  Array2 = %d.\n", nRemoved, Vec_IntSize(vArr1), Vec_IntSize(vArr2) );
    return -1;
}

void Gia_ManFindOneUnateInt( word * pOffSet, word * pOnSet, Vec_Ptr_t * vDivs, int nWords, Vec_Int_t * vUnateLits, Vec_Int_t * vNotUnateVars )
{
    word * pDiv; int i;
    Vec_IntClear( vUnateLits );
    Vec_IntClear( vNotUnateVars );
    Vec_PtrForEachEntryStart( word *, vDivs, pDiv, i, 2 )
        if ( !Abc_TtIntersectOne( pOffSet, 0, pDiv, 0, nWords ) )
            Vec_IntPush( vUnateLits, Abc_Var2Lit(i, 0) );
        else if ( !Abc_TtIntersectOne( pOffSet, 0, pDiv, 1, nWords ) )
            Vec_IntPush( vUnateLits, Abc_Var2Lit(i, 1) );
        else
            Vec_IntPush( vNotUnateVars, i );
}
int Gia_ManFindOneUnate( word * pSets[2], Vec_Ptr_t * vDivs, int nWords, Vec_Int_t * vUnateLits[2], Vec_Int_t * vNotUnateVars[2], int fVerbose )
{
    int n;
    for ( n = 0; n < 2; n++ )
    {
        Gia_ManFindOneUnateInt( pSets[n], pSets[!n], vDivs, nWords, vUnateLits[n], vNotUnateVars[n] );
        if ( fVerbose ) printf( "Found %d %d-unate divs.\n", Vec_IntSize(vUnateLits[n]), n );
    }
    return Gia_ManFindFirstCommonLit( vUnateLits[0], vUnateLits[1], fVerbose );
}

inline int Gia_ManDivCover( word * pOffSet, word * pOnSet, word * pDivA, int ComplA, word * pDivB, int ComplB, int nWords )
{
    //assert( !Abc_TtIntersectOne(pOffSet, 0, pDivA, ComplA, nWords) );
    //assert( !Abc_TtIntersectOne(pOffSet, 0, pDivB, ComplB, nWords) );
    return !Abc_TtIntersectTwo( pOnSet, 0, pDivA, !ComplA, pDivB, !ComplB, nWords );
}
int Gia_ManFindTwoUnateInt( word * pOffSet, word * pOnSet, Vec_Ptr_t * vDivs, int nWords, Vec_Int_t * vUnateLits )
{
    int i, k, iDiv0, iDiv1;
    Vec_IntForEachEntry( vUnateLits, iDiv1, i )
    Vec_IntForEachEntryStop( vUnateLits, iDiv0, k, i )
    {
        word * pDiv0 = (word *)Vec_PtrEntry(vDivs, Abc_Lit2Var(iDiv0));
        word * pDiv1 = (word *)Vec_PtrEntry(vDivs, Abc_Lit2Var(iDiv1));
        if ( Gia_ManDivCover(pOffSet, pOnSet, pDiv1, Abc_LitIsCompl(iDiv1), pDiv0, Abc_LitIsCompl(iDiv0), nWords) )
            return Abc_Var2Lit((Abc_LitNot(iDiv1) << 15) | Abc_LitNot(iDiv0), 1);
    }
    return -1;
}
int Gia_ManFindTwoUnate( word * pSets[2], Vec_Ptr_t * vDivs, int nWords, Vec_Int_t * vUnateLits[2], int fVerbose )
{
    int n, iLit;
    for ( n = 0; n < 2; n++ )
    {
        int nPairs = Vec_IntSize(vUnateLits[n])*(Vec_IntSize(vUnateLits[n])-1)/2;
        if ( fVerbose ) printf( "Trying %d pairs of %d-unate divs.\n", nPairs, n );
        iLit = Gia_ManFindTwoUnateInt( pSets[n], pSets[!n], vDivs, nWords, vUnateLits[n] );
        if ( iLit >= 0 )
            return Abc_LitNotCond(iLit, n);
    }
    return -1;
}

void Gia_ManFindXorInt( word * pOffSet, word * pOnSet, Vec_Int_t * vBinate, Vec_Ptr_t * vDivs, int nWords, Vec_Int_t * vUnatePairs )
{
    int i, k, iDiv0, iDiv1;
    Vec_IntForEachEntry( vBinate, iDiv1, i )
    Vec_IntForEachEntryStop( vBinate, iDiv0, k, i )
    {
        word * pDiv0 = (word *)Vec_PtrEntry(vDivs, iDiv0);
        word * pDiv1 = (word *)Vec_PtrEntry(vDivs, iDiv1);
        if ( !Abc_TtIntersectXor( pOffSet, 0, pDiv0, pDiv1, 0, nWords ) )
            Vec_IntPush( vUnatePairs, Abc_Var2Lit((Abc_Var2Lit(iDiv0, 0) << 15) | Abc_Var2Lit(iDiv1, 0), 0) );
        else if ( !Abc_TtIntersectXor( pOffSet, 0, pDiv0, pDiv1, 1, nWords ) )
            Vec_IntPush( vUnatePairs, Abc_Var2Lit((Abc_Var2Lit(iDiv0, 0) << 15) | Abc_Var2Lit(iDiv1, 0), 1) );
    }
}
int Gia_ManFindXor( word * pSets[2], Vec_Ptr_t * vDivs, int nWords, Vec_Int_t * vBinateVars, Vec_Int_t * vUnatePairs[2], int fVerbose )
{
    int n;
    for ( n = 0; n < 2; n++ )
    {
        Vec_IntClear( vUnatePairs[n] );
        Gia_ManFindXorInt( pSets[n], pSets[!n], vBinateVars, vDivs, nWords, vUnatePairs[n] );
        if ( fVerbose ) printf( "Found %d %d-unate XOR divs.\n", Vec_IntSize(vUnatePairs[n]), n );
    }
    return Gia_ManFindFirstCommonLit( vUnatePairs[0], vUnatePairs[1], fVerbose );
}

void Gia_ManFindUnatePairsInt( word * pOffSet, word * pOnSet, Vec_Int_t * vBinate, Vec_Ptr_t * vDivs, int nWords, Vec_Int_t * vUnatePairs )
{
    int n, i, k, iDiv0, iDiv1;
    Vec_IntForEachEntry( vBinate, iDiv1, i )
    Vec_IntForEachEntryStop( vBinate, iDiv0, k, i )
    {
        word * pDiv0 = (word *)Vec_PtrEntry(vDivs, iDiv0);
        word * pDiv1 = (word *)Vec_PtrEntry(vDivs, iDiv1);
        for ( n = 0; n < 4; n++ )
        {
            int iLit0 = Abc_Var2Lit( iDiv0, n&1 );
            int iLit1 = Abc_Var2Lit( iDiv1, n>>1 );
            if ( !Abc_TtIntersectTwo( pOffSet, 0, pDiv1, n>>1, pDiv0, n&1, nWords ) )
                Vec_IntPush( vUnatePairs, Abc_Var2Lit((iLit1 << 15) | iLit0, 0) );
        }
    }
}
void Gia_ManFindUnatePairs( word * pSets[2], Vec_Ptr_t * vDivs, int nWords, Vec_Int_t * vBinateVars, Vec_Int_t * vUnatePairs[2], int fVerbose )
{
    int n, RetValue;
    for ( n = 0; n < 2; n++ )
    {
        int nBefore = Vec_IntSize(vUnatePairs[n]);
        Gia_ManFindUnatePairsInt( pSets[n], pSets[!n], vBinateVars, vDivs, nWords, vUnatePairs[n] );
        if ( fVerbose ) printf( "Found %d %d-unate pair divs.\n", Vec_IntSize(vUnatePairs[n])-nBefore, n );
    }
    RetValue = Gia_ManFindFirstCommonLit( vUnatePairs[0], vUnatePairs[1], fVerbose );
    assert( RetValue == -1 );
}

void Gia_ManDeriveDivPair( int iDiv, Vec_Ptr_t * vDivs, int nWords, word * pRes )
{
    int fComp = Abc_LitIsCompl(iDiv);
    int iDiv0 = Abc_Lit2Var(iDiv) & 0x7FFF;
    int iDiv1 = Abc_Lit2Var(iDiv) >> 15;
    word * pDiv0 = (word *)Vec_PtrEntry(vDivs, Abc_Lit2Var(iDiv0));
    word * pDiv1 = (word *)Vec_PtrEntry(vDivs, Abc_Lit2Var(iDiv1));
    if ( iDiv0 < iDiv1 )
    {
        assert( !fComp );
        Abc_TtAndCompl( pRes, pDiv0, Abc_LitIsCompl(iDiv0), pDiv1, Abc_LitIsCompl(iDiv1), nWords );
    }
    else 
    {
        assert( !Abc_LitIsCompl(iDiv0) );
        assert( !Abc_LitIsCompl(iDiv1) );
        Abc_TtXor( pRes, pDiv0, pDiv1, nWords, 0 );
    }
}
int Gia_ManFindDivGateInt( word * pOffSet, word * pOnSet, Vec_Ptr_t * vDivs, int nWords, Vec_Int_t * vUnateLits, Vec_Int_t * vUnatePairs, word * pDivTemp )
{
    int i, k, iDiv0, iDiv1;
    int Limit1 = Abc_MinInt( Vec_IntSize(vUnateLits), 1000 );
    int Limit2 = Abc_MinInt( Vec_IntSize(vUnatePairs), 1000 );
    Vec_IntForEachEntryStop( vUnateLits, iDiv0, i, Limit1 )
    Vec_IntForEachEntryStop( vUnatePairs, iDiv1, k, Limit2 )
    {
        word * pDiv0 = (word *)Vec_PtrEntry(vDivs, Abc_Lit2Var(iDiv0));
        int fComp1   = Abc_LitIsCompl(iDiv1);
        Gia_ManDeriveDivPair( iDiv1, vDivs, nWords, pDivTemp );
        if ( Gia_ManDivCover(pOffSet, pOnSet, pDiv0, Abc_LitIsCompl(iDiv0), pDivTemp, fComp1, nWords) )
            return Abc_Var2Lit((Abc_Var2Lit(k, 1) << 15) | Abc_LitNot(iDiv0), 1);
    }
    return -1;
}
int Gia_ManFindDivGate( word * pSets[2], Vec_Ptr_t * vDivs, int nWords, Vec_Int_t * vUnateLits[2], Vec_Int_t * vUnatePairs[2], word * pDivTemp )
{
    int n, iLit;
    for ( n = 0; n < 2; n++ )
    {
        iLit = Gia_ManFindDivGateInt( pSets[n], pSets[!n], vDivs, nWords, vUnateLits[n], vUnatePairs[n], pDivTemp );
        if ( iLit >= 0 )
            return Abc_LitNotCond( iLit, n );
    }
    return -1;
}

int Gia_ManFindGateGateInt( word * pOffSet, word * pOnSet, Vec_Ptr_t * vDivs, int nWords, Vec_Int_t * vUnatePairs, word * pDivTempA, word * pDivTempB )
{
    int i, k, iDiv0, iDiv1;
    int Limit2 = Abc_MinInt( Vec_IntSize(vUnatePairs), 1000 );
    Vec_IntForEachEntryStop( vUnatePairs, iDiv1, i, Limit2 )
    {
        int fCompB = Abc_LitIsCompl(iDiv1);
        Gia_ManDeriveDivPair( iDiv1, vDivs, nWords, pDivTempB );
        Vec_IntForEachEntryStop( vUnatePairs, iDiv0, k, i )
        {
            int fCompA = Abc_LitIsCompl(iDiv0);
            Gia_ManDeriveDivPair( iDiv0, vDivs, nWords, pDivTempA );
            if ( Gia_ManDivCover(pOffSet, pOnSet, pDivTempA, fCompA, pDivTempB, fCompB, nWords) )
                return Abc_Var2Lit((Abc_Var2Lit(i, 1) << 15) | Abc_Var2Lit(k, 1), 1);
        }
    }
    return -1;
}
int Gia_ManFindGateGate( word * pSets[2], Vec_Ptr_t * vDivs, int nWords, Vec_Int_t * vUnatePairs[2], word * pDivTempA, word * pDivTempB )
{
    int n, iLit;
    for ( n = 0; n < 2; n++ )
    {
        iLit = Gia_ManFindGateGateInt( pSets[n], pSets[!n], vDivs, nWords, vUnatePairs[n], pDivTempA, pDivTempB );
        if ( iLit >= 0 )
            return Abc_LitNotCond( iLit, n );
    }
    return -1;
}

void Gia_ManComputeLitWeightsInt( word * pOffSet, word * pOnSet, Vec_Ptr_t * vDivs, int nWords, Vec_Int_t * vUnateLits, Vec_Int_t * vUnateLitsW )
{
    int i, iLit;
    Vec_IntClear( vUnateLitsW );
    Vec_IntForEachEntry( vUnateLits, iLit, i )
    {
        word * pDiv = (word *)Vec_PtrEntry( vDivs, Abc_Lit2Var(iLit) );
        //assert( !Abc_TtIntersectOne( pOffSet, 0, pDiv, Abc_LitIsCompl(iLit), nWords ) );
        Vec_IntPush( vUnateLitsW, -Abc_TtCountOnesVecMask(pDiv, pOnSet, nWords, Abc_LitIsCompl(iLit)) );
    }
}
void Gia_ManComputeLitWeights( word * pSets[2], Vec_Ptr_t * vDivs, int nWords, Vec_Int_t * vUnateLits[2], Vec_Int_t * vUnateLitsW[2], int TopW[2], int fVerbose )
{
    int n, TopNum = 5;
    TopW[0] = TopW[1] = 0;
    for ( n = 0; n < 2; n++ )
        Gia_ManComputeLitWeightsInt( pSets[n], pSets[!n], vDivs, nWords, vUnateLits[n], vUnateLitsW[n] );
    for ( n = 0; n < 2; n++ )
        if ( Vec_IntSize(vUnateLitsW[n]) )
        {
            int i, * pPerm = Abc_MergeSortCost( Vec_IntArray(vUnateLitsW[n]), Vec_IntSize(vUnateLitsW[n]) );
            TopW[n] = -Vec_IntEntry(vUnateLitsW[n], pPerm[0]);
            if ( fVerbose ) 
            {
                printf( "Top %d %d-unate divs:\n", TopNum, n );
                for ( i = 0; i < TopNum && i < Vec_IntSize(vUnateLits[n]); i++ )
                {
                    printf( "%5d : ", i );
                    printf( "Lit = %5d  ",  Vec_IntEntry(vUnateLits[n], pPerm[i]) );
                    printf( "Cost = %5d\n", -Vec_IntEntry(vUnateLitsW[n], pPerm[i]) );
                }
            }
            for ( i = 0; i < Vec_IntSize(vUnateLits[n]); i++ )
                pPerm[i] = Vec_IntEntry(vUnateLits[n], pPerm[i]);
            for ( i = 0; i < Vec_IntSize(vUnateLits[n]); i++ )
                vUnateLits[n]->pArray[i] = pPerm[i];
            ABC_FREE( pPerm );
        }
}

void Gia_ManComputePairWeightsInt( word * pOffSet, word * pOnSet, Vec_Ptr_t * vDivs, int nWords, Vec_Int_t * vUnatePairs, Vec_Int_t * vUnatePairsW )
{
    int i, iPair;
    Vec_IntClear( vUnatePairsW );
    Vec_IntForEachEntry( vUnatePairs, iPair, i )
    {
        int fComp = Abc_LitIsCompl(iPair);
        int iLit0 = Abc_Lit2Var(iPair) & 0x7FFF;
        int iLit1 = Abc_Lit2Var(iPair) >> 15;
        word * pDiv0 = (word *)Vec_PtrEntry( vDivs, Abc_Lit2Var(iLit0) );
        word * pDiv1 = (word *)Vec_PtrEntry( vDivs, Abc_Lit2Var(iLit1) );
        if ( iLit0 < iLit1 )
        {
            assert( !fComp );
            //assert( !Abc_TtIntersectTwo( pOffSet, 0, pDiv0, Abc_LitIsCompl(iLit0), pDiv1, Abc_LitIsCompl(iLit1), nWords ) );
            Vec_IntPush( vUnatePairsW, -Abc_TtCountOnesVecMask2(pDiv0, pDiv1, Abc_LitIsCompl(iLit0), Abc_LitIsCompl(iLit1), pOnSet, nWords) );
        }
        else
        {
            assert( !Abc_LitIsCompl(iLit0) );
            assert( !Abc_LitIsCompl(iLit1) );
            //assert( !Abc_TtIntersectXor( pOffSet, 0, pDiv0, pDiv1, fComp, nWords ) );
            Vec_IntPush( vUnatePairsW, -Abc_TtCountOnesVecXorMask(pDiv0, pDiv1, fComp, pOnSet, nWords) );
        }
    }
}
void Gia_ManComputePairWeights( word * pSets[2], Vec_Ptr_t * vDivs, int nWords, Vec_Int_t * vUnatePairs[2], Vec_Int_t * vUnatePairsW[2], int TopW[2], int fVerbose )
{
    int n, TopNum = 5;
    TopW[0] = TopW[1] = 0;
    for ( n = 0; n < 2; n++ )
        Gia_ManComputePairWeightsInt( pSets[n], pSets[!n], vDivs, nWords, vUnatePairs[n], vUnatePairsW[n] );
    for ( n = 0; n < 2; n++ )
        if ( Vec_IntSize(vUnatePairsW[n]) )
        {
            int i, * pPerm = Abc_MergeSortCost( Vec_IntArray(vUnatePairsW[n]), Vec_IntSize(vUnatePairsW[n]) );
            TopW[n] = -Vec_IntEntry(vUnatePairsW[n], pPerm[0]);
            if ( fVerbose ) 
            {
                printf( "Top %d %d-unate pairs:\n", TopNum, n );
                for ( i = 0; i < TopNum && i < Vec_IntSize(vUnatePairs[n]); i++ )
                {
                    int Pair = Vec_IntEntry(vUnatePairs[n], pPerm[i]);
                    int Div0 = Abc_Lit2Var(Pair) & 0x7FFF;
                    int Div1 = Abc_Lit2Var(Pair) >> 15;
                    printf( "%5d : ", i );
                    printf( "Compl = %5d  ", Abc_LitIsCompl(Pair) );
                    printf( "Type = %s  ",   Div0 < Div1 ? "and" : "xor" );
                    printf( "Div0 = %5d  ",  Div0 );
                    printf( "Div1 = %5d  ",  Div1 );
                    printf( "Cost = %5d\n",  -Vec_IntEntry(vUnatePairsW[n], pPerm[i]) );
                }
            }
            for ( i = 0; i < Vec_IntSize(vUnatePairs[n]); i++ )
                pPerm[i] = Vec_IntEntry(vUnatePairs[n], pPerm[i]);
            for ( i = 0; i < Vec_IntSize(vUnatePairs[n]); i++ )
                vUnatePairs[n]->pArray[i] = pPerm[i];
            ABC_FREE( pPerm );
        }
}

/**Function*************************************************************

  Synopsis    [Perform resubstitution.]

  Description []
               
  SideEffects []

  SeeAlso     []

***********************************************************************/
int Gia_ManResubPerform_rec( Gia_ResbMan_t * p, int nLimit )
{
    int TopOneW[2] = {0}, TopTwoW[2] = {0}, Max1, Max2, iResLit, nVars = Vec_PtrSize(p->vDivs);
    if ( p->fVerbose )
    {
        printf( "\nCalling decomposition for ISF: " ); 
        printf( "OFF = %5d (%6.2f %%)  ", Abc_TtCountOnesVec(p->pSets[0], p->nWords), 100.0*Abc_TtCountOnesVec(p->pSets[0], p->nWords)/(64*p->nWords) );
        printf( "ON = %5d (%6.2f %%)\n",  Abc_TtCountOnesVec(p->pSets[1], p->nWords), 100.0*Abc_TtCountOnesVec(p->pSets[1], p->nWords)/(64*p->nWords) );
    }
    if ( Abc_TtIsConst0( p->pSets[1], p->nWords ) )
        return 0;
    if ( Abc_TtIsConst0( p->pSets[0], p->nWords ) )
        return 1;
    iResLit = Gia_ManFindOneUnate( p->pSets, p->vDivs, p->nWords, p->vUnateLits, p->vNotUnateVars, p->fVerbose );
    if ( iResLit >= 0 ) // buffer
        return iResLit;
    if ( nLimit == 0 )
        return -1;
    iResLit = Gia_ManFindTwoUnate( p->pSets, p->vDivs, p->nWords, p->vUnateLits, p->fVerbose );
    if ( iResLit >= 0 ) // and
    {
        int iNode = nVars + Vec_IntSize(p->vGates)/2;
        int fComp = Abc_LitIsCompl(iResLit);
        int iDiv0 = Abc_Lit2Var(iResLit) & 0x7FFF;
        int iDiv1 = Abc_Lit2Var(iResLit) >> 15;
        assert( iDiv0 < iDiv1 );
        Vec_IntPushTwo( p->vGates, iDiv0, iDiv1 );
        return Abc_Var2Lit( iNode, fComp );
    }
    Vec_IntTwoFindCommon( p->vNotUnateVars[0], p->vNotUnateVars[1], p->vBinateVars );
    if ( Vec_IntSize(p->vBinateVars) > 1000 )
    {
        printf( "Reducing binate divs from %d to 1000.\n", Vec_IntSize(p->vBinateVars) );
        Vec_IntShrink( p->vBinateVars, 1000 );
    }
    iResLit = Gia_ManFindXor( p->pSets, p->vDivs, p->nWords, p->vBinateVars, p->vUnatePairs, p->fVerbose );
    if ( iResLit >= 0 ) // xor
    {
        int iNode = nVars + Vec_IntSize(p->vGates)/2;
        int fComp = Abc_LitIsCompl(iResLit);
        int iDiv0 = Abc_Lit2Var(iResLit) & 0x7FFF;
        int iDiv1 = Abc_Lit2Var(iResLit) >> 15;
        assert( !Abc_LitIsCompl(iDiv0) );
        assert( !Abc_LitIsCompl(iDiv1) );
        assert( iDiv0 > iDiv1 );
        Vec_IntPushTwo( p->vGates, iDiv0, iDiv1 );
        return Abc_Var2Lit( iNode, fComp );
    }
    if ( nLimit == 1 )
        return -1;
    Gia_ManFindUnatePairs( p->pSets, p->vDivs, p->nWords, p->vBinateVars, p->vUnatePairs, p->fVerbose );
    iResLit = Gia_ManFindDivGate( p->pSets, p->vDivs, p->nWords, p->vUnateLits, p->vUnatePairs, p->pDivA );
    if ( iResLit >= 0 ) // and(div,pair)
    {
        int iNode  = nVars + Vec_IntSize(p->vGates)/2;

        int fComp  = Abc_LitIsCompl(iResLit);
        int iDiv0  = Abc_Lit2Var(iResLit) & 0x7FFF; // div
        int iDiv1  = Abc_Lit2Var(iResLit) >> 15;    // pair

        int Div1   = Vec_IntEntry( p->vUnatePairs[!fComp], Abc_Lit2Var(iDiv1) );
        int fComp1 = Abc_LitIsCompl(Div1) ^ Abc_LitIsCompl(iDiv1);
        int iDiv10 = Abc_Lit2Var(Div1) & 0x7FFF;
        int iDiv11 = Abc_Lit2Var(Div1) >> 15;   

        Vec_IntPushTwo( p->vGates, iDiv10, iDiv11 );
        Vec_IntPushTwo( p->vGates, iDiv0, Abc_Var2Lit(iNode, fComp1) );
        return Abc_Var2Lit( iNode+1, fComp );
    }
    if ( nLimit == 2 )
        return -1;
    iResLit = Gia_ManFindGateGate( p->pSets, p->vDivs, p->nWords, p->vUnatePairs, p->pDivA, p->pDivB );
    if ( iResLit >= 0 ) // and(pair,pair)
    {
        int iNode  = nVars + Vec_IntSize(p->vGates)/2;

        int fComp  = Abc_LitIsCompl(iResLit);
        int iDiv0  = Abc_Lit2Var(iResLit) & 0x7FFF; // pair
        int iDiv1  = Abc_Lit2Var(iResLit) >> 15;    // pair

        int Div0   = Vec_IntEntry( p->vUnatePairs[!fComp], Abc_Lit2Var(iDiv0) );
        int fComp0 = Abc_LitIsCompl(Div0) ^ Abc_LitIsCompl(iDiv0);
        int iDiv00 = Abc_Lit2Var(Div0) & 0x7FFF;
        int iDiv01 = Abc_Lit2Var(Div0) >> 15;   
        
        int Div1   = Vec_IntEntry( p->vUnatePairs[!fComp], Abc_Lit2Var(iDiv1) );
        int fComp1 = Abc_LitIsCompl(Div1) ^ Abc_LitIsCompl(iDiv1);
        int iDiv10 = Abc_Lit2Var(Div1) & 0x7FFF;
        int iDiv11 = Abc_Lit2Var(Div1) >> 15;   
        
        Vec_IntPushTwo( p->vGates, iDiv00, iDiv01 );
        Vec_IntPushTwo( p->vGates, iDiv10, iDiv11 );
        Vec_IntPushTwo( p->vGates, Abc_Var2Lit(iNode, fComp0), Abc_Var2Lit(iNode+1, fComp1) );
        return Abc_Var2Lit( iNode+2, fComp );
    }
    if ( nLimit == 3 )
        return -1;
    if ( Vec_IntSize(p->vUnateLits[0]) + Vec_IntSize(p->vUnateLits[1]) + Vec_IntSize(p->vUnatePairs[0]) + Vec_IntSize(p->vUnatePairs[1]) == 0 )
        return -1;
    Gia_ManComputeLitWeights( p->pSets, p->vDivs, p->nWords, p->vUnateLits, p->vUnateLitsW, TopOneW, p->fVerbose );
    Gia_ManComputePairWeights( p->pSets, p->vDivs, p->nWords, p->vUnatePairs, p->vUnatePairsW, TopTwoW, p->fVerbose );
    Max1 = Abc_MaxInt(TopOneW[0], TopOneW[1]);
    Max2 = Abc_MaxInt(TopTwoW[0], TopTwoW[1]);
    if ( Abc_MaxInt(Max1, Max2) == 0 )
        return -1;
    if ( Max1 > Max2/2 )
    {
        if ( Max1 == TopOneW[0] || Max1 == TopOneW[1] )
        {
            int fUseOr  = Max1 == TopOneW[0];
            int iDiv    = Vec_IntEntry( p->vUnateLits[!fUseOr], 0 );
            int fComp   = Abc_LitIsCompl(iDiv);
            word * pDiv = (word *)Vec_PtrEntry( p->vDivs, Abc_Lit2Var(iDiv) );
            Abc_TtAndSharp( p->pSets[fUseOr], p->pSets[fUseOr], pDiv, p->nWords, !fComp );
            iResLit = Gia_ManResubPerform_rec( p, nLimit-1 );
            if ( iResLit >= 0 ) 
            {
                int iNode = nVars + Vec_IntSize(p->vGates)/2;
                Vec_IntPushTwo( p->vGates, Abc_LitNot(iDiv), Abc_LitNotCond(iResLit, fUseOr) );
                return Abc_Var2Lit( iNode, fUseOr );
            }
        }
        if ( Max2 == 0 )
            return -1;
        if ( Max2 == TopTwoW[0] || Max2 == TopTwoW[1] )
        {
            int fUseOr  = Max2 == TopTwoW[0];
            int iDiv    = Vec_IntEntry( p->vUnatePairs[!fUseOr], 0 );
            int fComp   = Abc_LitIsCompl(iDiv);
            Gia_ManDeriveDivPair( iDiv, p->vDivs, p->nWords, p->pDivA );
            Abc_TtAndSharp( p->pSets[fUseOr], p->pSets[fUseOr], p->pDivA, p->nWords, !fComp );
            iResLit = Gia_ManResubPerform_rec( p, nLimit-2 );
            if ( iResLit >= 0 ) 
            {
                int iNode = nVars + Vec_IntSize(p->vGates)/2;
                int iDiv0 = Abc_Lit2Var(iDiv) & 0x7FFF;   
                int iDiv1 = Abc_Lit2Var(iDiv) >> 15;      
                Vec_IntPushTwo( p->vGates, iDiv0, iDiv1 );
                Vec_IntPushTwo( p->vGates, Abc_LitNotCond(iResLit, fUseOr), Abc_Var2Lit(iNode, !fComp) );
                return Abc_Var2Lit( iNode+1, fUseOr );
            }
        }
    }
    else
    {
        if ( Max2 == TopTwoW[0] || Max2 == TopTwoW[1] )
        {
            int fUseOr  = Max2 == TopTwoW[0];
            int iDiv    = Vec_IntEntry( p->vUnatePairs[!fUseOr], 0 );
            int fComp   = Abc_LitIsCompl(iDiv);
            Gia_ManDeriveDivPair( iDiv, p->vDivs, p->nWords, p->pDivA );
            Abc_TtAndSharp( p->pSets[fUseOr], p->pSets[fUseOr], p->pDivA, p->nWords, !fComp );
            iResLit = Gia_ManResubPerform_rec( p, nLimit-2 );
            if ( iResLit >= 0 ) 
            {
                int iNode = nVars + Vec_IntSize(p->vGates)/2;
                int iDiv0 = Abc_Lit2Var(iDiv) & 0x7FFF;   
                int iDiv1 = Abc_Lit2Var(iDiv) >> 15;      
                Vec_IntPushTwo( p->vGates, iDiv0, iDiv1 );
                Vec_IntPushTwo( p->vGates, Abc_LitNotCond(iResLit, fUseOr), Abc_Var2Lit(iNode, !fComp) );
                return Abc_Var2Lit( iNode+1, fUseOr );
            }
        }
        if ( Max1 == 0 )
            return -1;
        if ( Max1 == TopOneW[0] || Max1 == TopOneW[1] )
        {
            int fUseOr  = Max1 == TopOneW[0];
            int iDiv    = Vec_IntEntry( p->vUnateLits[!fUseOr], 0 );
            int fComp   = Abc_LitIsCompl(iDiv);
            word * pDiv = (word *)Vec_PtrEntry( p->vDivs, Abc_Lit2Var(iDiv) );
            Abc_TtAndSharp( p->pSets[fUseOr], p->pSets[fUseOr], pDiv, p->nWords, !fComp );
            iResLit = Gia_ManResubPerform_rec( p, nLimit-1 );
            if ( iResLit >= 0 ) 
            {
                int iNode = nVars + Vec_IntSize(p->vGates)/2;
                Vec_IntPushTwo( p->vGates, Abc_LitNot(iDiv), Abc_LitNotCond(iResLit, fUseOr) );
                return Abc_Var2Lit( iNode, fUseOr );
            }
        }
    }
    return -1;
}
void Gia_ManResubPerform( Gia_ResbMan_t * p, Vec_Ptr_t * vDivs, int nWords, int nLimit, int fDebug, int fVerbose )
{
    int Res;
    Gia_ResbInit( p, vDivs, nWords, nLimit, fDebug, fVerbose );
    Res = Gia_ManResubPerform_rec( p, nLimit );
    if ( Res >= 0 )
        Vec_IntPush( p->vGates, Res );
}

/**Function*************************************************************

  Synopsis    [Top level.]

  Description []
               
  SideEffects []

  SeeAlso     []

***********************************************************************/
static Gia_ResbMan_t * s_pResbMan = NULL;

void Abc_ResubPrepareManager( int nWords )
{
    if ( s_pResbMan != nullptr )
    {
      Gia_ResbFree( s_pResbMan );
      s_pResbMan = nullptr;
    }
    if ( nWords > 0 )
    {
      s_pResbMan = Gia_ResbAlloc( nWords );
    }
}

int Abc_ResubComputeFunction( void ** ppDivs, int nDivs, int nWords, int nLimit, int fDebug, int fVerbose, int ** ppArray )
{
    Vec_Ptr_t Divs = { nDivs, nDivs, ppDivs };
    assert( s_pResbMan != NULL ); // first call Abc_ResubPrepareManager()
    Gia_ManResubPerform( s_pResbMan, &Divs, nWords, nLimit, fDebug, fVerbose );
    if ( fVerbose )
    {
        Gia_ManResubPrint( s_pResbMan->vGates, nDivs );
        printf( "\n" );
    }
    if ( fDebug )
    {
        if ( !Gia_ManResubVerify(s_pResbMan) )
        {
            Gia_ManResubPrint( s_pResbMan->vGates, nDivs );
            printf( "Verification FAILED.\n" );
        }
        //else
        //    printf( "Verification succeeded.\n" );
    }
    *ppArray = Vec_IntArray(s_pResbMan->vGates);
    assert( Vec_IntSize(s_pResbMan->vGates)/2 <= nLimit );
    return Vec_IntSize(s_pResbMan->vGates);
}

#if 0
void Abc_ResubDumpProblem( char * pFileName, void ** ppDivs, int nDivs, int nWords )
{
    Vec_Wrd_t * vSims = Vec_WrdAlloc( nDivs * nWords );
    word ** pDivs = (word **)ppDivs;
    int d, w;
    for ( d = 0; d < nDivs;  d++ )
    for ( w = 0; w < nWords; w++ )
        Vec_WrdPush( vSims, pDivs[d][w] );
    Gia_ManSimPatWrite( pFileName, vSims, nWords );
    Vec_WrdFree( vSims );
}
#endif

} /* namespace abcresub */
