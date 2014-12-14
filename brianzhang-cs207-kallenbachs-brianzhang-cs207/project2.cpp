#include <fstream>

#include "CS207/SDLViewer.hpp"
#include "CS207/Util.hpp"
#include "CS207/Color.hpp"

#include "Graph.hpp"
#include "Point.hpp"
#include "Solid.hpp"
#include "../jb_parallel.hpp"


// Gravity in meters/sec^2
static constexpr double grav = 9.81;

int main() {
  Interactables::Ball b(Point(0, 0, 5), 1, Point(0, 0, 0), 1);
  Interactables::Floor f(0, 0.9);

  Interactables its;

  // CS207::SDLViewer viewer;
  // viewer.launch();
  // b.draw(viewer);
  // f.draw(viewer);
  // viewer.center_view();

  double dt = 0.001;
  double t_start = 0.0;
  double t_end   = 3.0;

  { jb_parallel::Timer timer("");
  for (double t = t_start; t < t_end; t += dt) {
    b.step(grav, t, dt);
    its.interact(b, f);
    // viewer.clear();
    // b.draw(viewer);
    // f.draw(viewer);
    // viewer.set_label(t);
  }
  }

  return 0;
}
