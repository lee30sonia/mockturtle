/*!
  \file tt.hpp
  \brief Extracted from ABC
         https://github.com/berkeley-abc/abc

  \author Alan Mishchenko (UC Berkeley)
*/

#pragma once

namespace abcresub
{

inline int Abc_TtCountOnes( word x )
{
    x = x - ((x >> 1) & ABC_CONST(0x5555555555555555));   
    x = (x & ABC_CONST(0x3333333333333333)) + ((x >> 2) & ABC_CONST(0x3333333333333333));    
    x = (x + (x >> 4)) & ABC_CONST(0x0F0F0F0F0F0F0F0F);    
    x = x + (x >> 8);
    x = x + (x >> 16);
    x = x + (x >> 32); 
    return (int)(x & 0xFF);
}

inline int Abc_TtCountOnesVec( word * x, int nWords )
{
    int w, Count = 0;
    for ( w = 0; w < nWords; w++ )
        Count += Abc_TtCountOnes( x[w] );
    return Count;
}

inline int Abc_TtCountOnesVecMask( word * x, word * pMask, int nWords, int fCompl )
{
    int w, Count = 0;
    if ( fCompl )
        for ( w = 0; w < nWords; w++ )
            Count += Abc_TtCountOnes( pMask[w] & ~x[w] );
    else
        for ( w = 0; w < nWords; w++ )
            Count += Abc_TtCountOnes( pMask[w] & x[w] );
    return Count;
}

inline int Abc_TtCountOnesVecMask2( word * x0, word * x1, int fComp0, int fComp1, word * pMask, int nWords )
{
    int w, Count = 0;
    if ( !fComp0 && !fComp1 )
        for ( w = 0; w < nWords; w++ )
            Count += Abc_TtCountOnes( pMask[w] &  x0[w] &  x1[w] );
    else if (  fComp0 && !fComp1 )
        for ( w = 0; w < nWords; w++ )
            Count += Abc_TtCountOnes( pMask[w] & ~x0[w] &  x1[w] );
    else if ( !fComp0 &&  fComp1 )
        for ( w = 0; w < nWords; w++ )
            Count += Abc_TtCountOnes( pMask[w] &  x0[w] & ~x1[w] );
    else 
        for ( w = 0; w < nWords; w++ )
            Count += Abc_TtCountOnes( pMask[w] & ~x0[w] & ~x1[w] );
    return Count;
}

inline int Abc_TtCountOnesVecXorMask( word * x, word * y, int fCompl, word * pMask, int nWords )
{
    int w, Count = 0;
    if ( fCompl )
        for ( w = 0; w < nWords; w++ )
            Count += Abc_TtCountOnes( pMask[w] & (x[w] ^ ~y[w]) );
    else
        for ( w = 0; w < nWords; w++ )
            Count += Abc_TtCountOnes( pMask[w] & (x[w] ^ y[w]) );
    return Count;
}

inline void Abc_TtCopy( word * pOut, word * pIn, int nWords, int fCompl )
{
    int w;
    if ( fCompl )
        for ( w = 0; w < nWords; w++ )
            pOut[w] = ~pIn[w];
    else
        for ( w = 0; w < nWords; w++ )
            pOut[w] = pIn[w];
}

inline int Abc_TtIsConst0( word * pIn1, int nWords )
{
    int w;
    for ( w = 0; w < nWords; w++ )
        if ( pIn1[w] )
            return 0;
    return 1;
}

inline void Abc_TtClear( word * pOut, int nWords )
{
    int w;
    for ( w = 0; w < nWords; w++ )
        pOut[w] = 0;
}

inline void Abc_TtFill( word * pOut, int nWords )
{
    int w;
    for ( w = 0; w < nWords; w++ )
        pOut[w] = ~(word)0;
}

inline void Abc_TtAndCompl( word * pOut, word * pIn1, int fCompl1, word * pIn2, int fCompl2, int nWords )
{
    int w;
    if ( fCompl1 )
    {
        if ( fCompl2 )
            for ( w = 0; w < nWords; w++ )
                pOut[w] = ~pIn1[w] & ~pIn2[w];
        else
            for ( w = 0; w < nWords; w++ )
                pOut[w] = ~pIn1[w] & pIn2[w];
    }
    else
    {
        if ( fCompl2 )
            for ( w = 0; w < nWords; w++ )
                pOut[w] = pIn1[w] & ~pIn2[w];
        else
            for ( w = 0; w < nWords; w++ )
                pOut[w] = pIn1[w] & pIn2[w];
    }
}

inline void Abc_TtAndSharp( word * pOut, word * pIn1, word * pIn2, int nWords, int fCompl )
{
    int w;
    if ( fCompl )
        for ( w = 0; w < nWords; w++ )
            pOut[w] = pIn1[w] & ~pIn2[w];
    else
        for ( w = 0; w < nWords; w++ )
            pOut[w] = pIn1[w] & pIn2[w];
}

inline void Abc_TtXor( word * pOut, word * pIn1, word * pIn2, int nWords, int fCompl )
{
    int w;
    if ( fCompl )
        for ( w = 0; w < nWords; w++ )
            pOut[w] = pIn1[w] ^ ~pIn2[w];
    else
        for ( w = 0; w < nWords; w++ )
            pOut[w] = pIn1[w] ^ pIn2[w];
}

inline int Abc_TtIntersectOne( word * pOut, int fComp, word * pIn, int fComp0, int nWords )
{
    int w;
    if ( fComp0 )
    {
        if ( fComp )
        {
            for ( w = 0; w < nWords; w++ )
                if ( ~pIn[w] & ~pOut[w] )
                    return 1;
        }
        else
        {
            for ( w = 0; w < nWords; w++ )
                if ( ~pIn[w] & pOut[w] )
                    return 1;
        }
    }
    else
    {
        if ( fComp )
        {
            for ( w = 0; w < nWords; w++ )
                if ( pIn[w] & ~pOut[w] )
                    return 1;
        }
        else
        {
            for ( w = 0; w < nWords; w++ )
                if ( pIn[w] & pOut[w] )
                    return 1;
        }
    }
    return 0;
}

inline int Abc_TtIntersectTwo( word * pOut, int fComp, word * pIn0, int fComp0, word * pIn1, int fComp1, int nWords )
{
    int w;
    if ( fComp0 && fComp1 )
    {
        if ( fComp )
        {
            for ( w = 0; w < nWords; w++ )
                if ( ~pIn0[w] & ~pIn1[w] & ~pOut[w] )
                    return 1;
        }
        else
        {
            for ( w = 0; w < nWords; w++ )
                if ( ~pIn0[w] & ~pIn1[w] & pOut[w] )
                    return 1;
        }
    }
    else if ( fComp0 )
    {
        if ( fComp )
        {
            for ( w = 0; w < nWords; w++ )
                if ( ~pIn0[w] & pIn1[w] & ~pOut[w] )
                    return 1;
        }
        else
        {
            for ( w = 0; w < nWords; w++ )
                if ( ~pIn0[w] & pIn1[w] & pOut[w] )
                    return 1;
        }
    }
    else if ( fComp1 )
    {
        if ( fComp )
        {
            for ( w = 0; w < nWords; w++ )
                if ( pIn0[w] & ~pIn1[w] & ~pOut[w] )
                    return 1;
        }
        else
        {
            for ( w = 0; w < nWords; w++ )
                if ( pIn0[w] & ~pIn1[w] & pOut[w] )
                    return 1;
        }
    }
    else 
    {
        if ( fComp )
        {
            for ( w = 0; w < nWords; w++ )
                if ( pIn0[w] & pIn1[w] & ~pOut[w] )
                    return 1;
        }
        else
        {
            for ( w = 0; w < nWords; w++ )
                if ( pIn0[w] & pIn1[w] & pOut[w] )
                    return 1;
        }
    }
    return 0;
}

inline int Abc_TtIntersectXor( word * pOut, int fComp, word * pIn0, word * pIn1, int fComp01, int nWords )
{
    int w;
    if ( fComp01 )
    {
        if ( fComp )
        {
            for ( w = 0; w < nWords; w++ )
                if ( ~(pIn0[w] ^ pIn1[w]) & ~pOut[w] )
                    return 1;
        }
        else 
        {
            for ( w = 0; w < nWords; w++ )
                if ( ~(pIn0[w] ^ pIn1[w]) & pOut[w] )
                    return 1;
        }
    }
    else
    {
        if ( fComp )
        {
            for ( w = 0; w < nWords; w++ )
                if ( (pIn0[w] ^ pIn1[w]) & ~pOut[w] )
                    return 1;
        }
        else 
        {
            for ( w = 0; w < nWords; w++ )
                if ( (pIn0[w] ^ pIn1[w]) & pOut[w] )
                    return 1;
        }
    }
    return 0;
}

} /* abcresub */
