#include <math.h>

#include <SDL/SDL.h>

#include "live_view.h"

using namespace LAMMPS_NS;

LiveView::LiveView(int in_nx, int in_ny, int in_sx, int in_sy, Error *error)
{
  nx_ = in_nx;
  ny_ = in_ny;

  sx_ = in_sx;
  sy_ = in_sy;

  if (SDL_Init(SDL_INIT_VIDEO))
    error->all("Could not initialize SDL.");
  screen_ = SDL_SetVideoMode(sx_, sy_, 0, 0);
}


LiveView::~LiveView()
{
  SDL_Quit();
}


template<typename T>
void update_live_view(SDL_Surface *screen, int nx, int ny, int sx, int sy,
		      double *map, T color0, T colorrep, T coloratt)
{
  for (int j = 0; j < sy; j++) {
    T *bufp = ((T*) screen->pixels) + j*screen->pitch/sizeof(T);
    for (int i = 0; i < sx; i++, bufp++) {
      double f = map[j*ny*nx/sy + i*nx/sx];
      if (f < 0.0) {
	*bufp = colorrep;
      }
      else if (f > 0.0) {
	*bufp = coloratt;
      }
      else {
	*bufp = color0;
      }
    }
  }
}


void LiveView::update(double *map)
{
  /* Allow to move window */
  SDL_Event event;
  SDL_PollEvent(&event);

  if (SDL_MUSTLOCK(screen_)) {
    if (SDL_LockSurface(screen_) < 0)
      return;
  }

  Uint32 black = SDL_MapRGB(screen_->format, 0, 0, 0);
  Uint32 white = SDL_MapRGB(screen_->format, 255, 255, 255);
  Uint32 orange = SDL_MapRGB(screen_->format, 255, 127, 0);

  switch (screen_->format->BytesPerPixel) {
  case 1:
    update_live_view<Uint8>(screen_, nx_, ny_, sx_, sy_, map,
			    black, orange, white);
    break;

  case 2:
    update_live_view<Uint16>(screen_, nx_, ny_, sx_, sy_, map,
			     black, orange, white);
    break;

  case 4:
    update_live_view<Uint32>(screen_, nx_, ny_, sx_, sy_, map,
			     black, orange, white);
    break;
  }

  if (SDL_MUSTLOCK(screen_)) {
    SDL_UnlockSurface(screen_);
  }

  SDL_Flip(screen_);
}

