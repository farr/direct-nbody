open OUnit

module Ic = Ics.Make(Leapfrog)
module An = Analysis.Make(Leapfrog)

let test_neighbors () = 
  let bs = Ic.make_plummer 1000 in 
  let b = bs.(101) and 
      t = An.bodies_to_neighbor_tree bs in 
  let nbrs = An.nearest_neighbors 6 b t and 
      nbrs' = An.nearest_neighbors_slow 6 bs b in 
    Array.fast_sort Pervasives.compare nbrs;
    Array.fast_sort Pervasives.compare nbrs';
    for i = 0 to 5 do 
      assert_equal nbrs.(i) nbrs'.(i)
    done

let test_density_radius () = 
  let bs = Ic.make_plummer 10000 in 
  let rhos = An.density_squared_estimators 6 bs in 
  let r = An.density_radius rhos bs in 
  let rexact = 189.0 /. 640.0 in
    Printf.fprintf stderr "r = %g\n" r;
    assert_bool "r too large" (r < 2.0*.rexact);
    assert_bool "r too small" (r > 0.5*.rexact)

let tests = "analysis.ml tests" >:::
  ["nearest neighbor tests" >:: test_neighbors;
   "density radius test" >:: test_density_radius]
