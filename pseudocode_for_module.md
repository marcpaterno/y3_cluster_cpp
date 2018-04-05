## Pseudocode for CosmoSIS module class template

This pseudocode makes several simplifications. Primary among them is the
reduction of the number of integrand terms from 12 to 2, with a concomitant
reduction in the number of class template parameters.

First, here is the pseudocode for the class template:

```c++

#include "gamma_t.hh"
... lots of other headers probably needed!

template <class TERM1, class TERM2>
class Y3_Clustering {
private:
  cubacpp::Cuhre integrator_; // the apparent champ of integrators

  // Constructor from a cosmosis::DataBlock.
  // See cosmosis/datablock/datablock.hh for the description
  // of the class DataBlock.
  explicit Y3_Clustering(cosmosis::DataBlock const& config)
    : cuhre_()
    , epsrel_(),
    , epsabs_()
    {
      config.get_val("???", "maxeval", cuhre_.maxeval); // set cuhre_.maxeval
      config.get_val("???", "rule", cuhre_.key); // set rule for CUHRE 
      config.get_val("???", "epsrel", epsrel_);
      config.get_val("???", "epsabs", epsabs_);
    }
  

  // This is the function that will be called on every sample.
  void execute(cosmosis::DataBlock& sample) const {
    // Construct the integrand object
    TERM1 t1(sample); // create the term, using info from the current sample
    TERM2 t2(sample);
    auto igrand = make_gamma_t_integrand(t1, t2);
    // ... and integrate ...
    auto res = integrator_.integrate(igrand, epsrel_, epsabs_);
    // now put the result of the integration, and any other information we care
    // to record, into the datablock...
    sample.put_val("???", "good_name_here", res.value);
    ...
  }
}
```

Assume the class template is in the file `y3_clustering.hh`. The user writes
as many headers as he wants, containing term implementations. Maybe we have
a bunch of "standard" ones in the repository. He also writes a single `.cc` file:

```c++

// File MyY3Module.cc

#include "y3_clustering.hh"
#include "term1.hh" // Use this implementation for term1
#include "term2.hh" // ... and this for term2

// Introduce a type alias ("typedef") that is the module class I want to use:

using MyY3Module = y3_clutsering<term1, term2>;

// Use the macro that generates the code to make this a CosmoSIS module:

DEFINE_COSMOSIS_MODULE(MyY3Module);
```

The user compiles `MyY3Module.cc`, links it with the right libraries, and
produces a dynamic library `MyY3Module.so`. This library is the CosmoSIS
module.
