// $Header: /nfs/dsi/cvs/galib/ga/GASimpleGA.h,v 1.1.1.1 1999/11/11 18:56:03 mbwall Exp $
/* ----------------------------------------------------------------------------
  gasimple.h
  mbwall 28jul94
  Copyright (c) 1995 Massachusetts Institute of Technology
                     all rights reserved

  Header file for the simple genetic algorithm class.
---------------------------------------------------------------------------- */
#ifndef _ga_gasimple_h_
#define _ga_gasimple_h_

#include <ga/GABaseGA.h>

class GASimpleGA : public GAGeneticAlgorithm {
public:
  GADefineIdentity("GASimpleGA", GAID::SimpleGA);

  static GAParameterList& registerDefaultParameters(GAParameterList&);

public:
  GASimpleGA(const GAGenome&);
  GASimpleGA(const GAPopulation&);
  GASimpleGA(const GASimpleGA&);
  GASimpleGA& operator=(const GASimpleGA&);
  virtual ~GASimpleGA();
  virtual void copy(const GAGeneticAlgorithm&);

  virtual void initialize(unsigned int seed=0);
  virtual void step();
  GASimpleGA & operator++() { step(); return *this; }

  virtual int setptr(const char* name, const void* value);
  virtual int get(const char* name, void* value) const;

  GABoolean elitist() const {return el;}
  GABoolean elitist(GABoolean flag)
    {params.set(gaNelitism, (int)flag); return el=flag;}

  virtual int minimaxi() const {return minmax;}
  virtual int minimaxi(int m);

  virtual const GAPopulation& population() const {return *pop;}
  virtual const GAPopulation& population(const GAPopulation&);
  virtual int populationSize() const {return pop->size();}
  virtual int populationSize(unsigned int n);
  virtual GAScalingScheme& scaling() const {return pop->scaling();}
  virtual GAScalingScheme& scaling(const GAScalingScheme & s)
    {oldPop->scaling(s); return GAGeneticAlgorithm::scaling(s);}
  virtual GASelectionScheme& selector() const {return pop->selector(); }
  virtual GASelectionScheme& selector(const GASelectionScheme& s)
    {oldPop->selector(s); return GAGeneticAlgorithm::selector(s);}
  virtual void objectiveFunction(GAGenome::Evaluator f);
  virtual void objectiveData(const GAEvalData& v);

protected:
  GAPopulation *oldPop;		// current and old populations
  GABoolean el;			// are we elitist?
};



#ifndef NO_STREAMS
inline ostream& operator<< (ostream& os, GASimpleGA & arg)
{ arg.write(os); return(os); }
inline istream& operator>> (istream& is, GASimpleGA & arg)
{ arg.read(is); return(is); }
#endif

#endif
