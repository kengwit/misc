#ifndef NZ_H
# define NZ_H

template <int NUM_COMPONENTS>
class NZ {

 public:

  NZ (FIELD<NUM_COMPONENTS> *f) : _f(f) { }
  size_t nzs (BFUN_DOFPARAM_PAIR_ID dpid);

 private:

  FIELD<NUM_COMPONENTS> *_f;

};

#endif
