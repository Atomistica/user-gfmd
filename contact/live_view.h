#ifndef __LIVE_VIEW
#define __LIVE_VIEW

#include "error.h"

namespace LAMMPS_NS {

class LiveView {
 public:
  LiveView(int, int, int, int, Error *);
  ~LiveView();

  void update(double *);

 protected:
  int nx_, ny_, sx_, sy_;

  class SDL_Surface *screen_;
};

}


#endif
