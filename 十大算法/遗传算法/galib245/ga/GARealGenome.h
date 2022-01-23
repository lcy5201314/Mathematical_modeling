// $Header: /nfs/dsi/cvs/galib/ga/GARealGenome.h,v 1.1.1.1 1999/11/11 18:56:03 mbwall Exp $
/* ----------------------------------------------------------------------------
  real.h
  mbwall 25feb95
  Copyright (c) 1995 Massachusetts Institute of Technology
                     all rights reserved

 DESCRIPTION:
  This header defines the specialization of the array genome of type float
for the real number genome.
---------------------------------------------------------------------------- */
#ifndef _ga_real_h_
#define _ga_real_h_

#include <ga/GAAllele.h>
#include <ga/GA1DArrayGenome.h>

typedef GAAlleleSet<float> GARealAlleleSet;
typedef GAAlleleSetArray<float> GARealAlleleSetArray;

typedef GA1DArrayAlleleGenome<float> GARealGenome;

inline void GARealUniformInitializer(GAGenome& g){
  GA1DArrayAlleleGenome<float>::UniformInitializer(g);
}
inline void GARealOrderedInitializer(GAGenome& g){
  GA1DArrayAlleleGenome<float>::OrderedInitializer(g);
}

inline int GARealUniformMutator(GAGenome& g, float pmut){
  return GA1DArrayAlleleGenome<float>::FlipMutator(g, pmut);
}
inline int GARealSwapMutator(GAGenome& g, float pmut){
  return GA1DArrayGenome<float>::SwapMutator(g, pmut);
}
int GARealGaussianMutator(GAGenome &, float);


inline int GARealUniformCrossover(const GAGenome& a, const GAGenome& b,
				  GAGenome* c, GAGenome* d) {
  return GA1DArrayGenome<float>::UniformCrossover(a,b,c,d);
}
inline int GARealEvenOddCrossover(const GAGenome& a, const GAGenome& b,
				  GAGenome* c, GAGenome* d) {
  return GA1DArrayGenome<float>::EvenOddCrossover(a,b,c,d);
}
inline int GARealOnePointCrossover(const GAGenome& a, const GAGenome& b,
				  GAGenome* c, GAGenome* d) {
  return GA1DArrayGenome<float>::OnePointCrossover(a,b,c,d);
}
inline int GARealTwoPointCrossover(const GAGenome& a, const GAGenome& b,
				  GAGenome* c, GAGenome* d) {
  return GA1DArrayGenome<float>::TwoPointCrossover(a,b,c,d);
}
inline int GARealPartialMatchCrossover(const GAGenome& a, const GAGenome& b,
				  GAGenome* c, GAGenome* d) {
  return GA1DArrayGenome<float>::PartialMatchCrossover(a,b,c,d);
}
inline int GARealOrderCrossover(const GAGenome& a, const GAGenome& b,
				  GAGenome* c, GAGenome* d) {
  return GA1DArrayGenome<float>::OrderCrossover(a,b,c,d);
}
inline int GARealCycleCrossover(const GAGenome& a, const GAGenome& b,
				  GAGenome* c, GAGenome* d) {
  return GA1DArrayGenome<float>::CycleCrossover(a,b,c,d);
}
int GARealArithmeticCrossover(const GAGenome& a, const GAGenome& b,
				  GAGenome* c, GAGenome* d);
int GARealBlendCrossover(const GAGenome& a, const GAGenome& b,
				  GAGenome* c, GAGenome* d);

#endif
