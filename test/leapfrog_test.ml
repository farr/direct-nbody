open OUnit
open Leapfrog

module Ic = Ics.Make(Leapfrog)
module E = Energy.Make(Leapfrog)

let test_const_dt_energy () = 
  let bs = Ic.rescale_to_standard_units (Ic.adjust_frame (Ic.make_hot_spherical 100)) in 
  let bs1 = advance_const_dt bs 1e-4 1 and 
      bs2 = advance_const_dt bs 5e-5 1 in 
  let de1 = abs_float (((E.energy bs1) +. 0.25) /. 0.25) and 
      de2 = abs_float (((E.energy bs2) +. 0.25) /. 0.25) in 
  let r = de1 /. de2 in 
    assert_bool (Printf.sprintf "energy scaling = %g (should be about 8)" r)
      ((r < (sqrt (8.0*.16.0))) && (r > (sqrt (4.0*.8.0))))

let test_leap_energy () = 
  let bs = Ic.rescale_to_standard_units (Ic.adjust_frame (Ic.make_hot_spherical 100)) in 
  let e1 = E.energy (advance bs 1.0 1e-3) and
      e2 = E.energy (advance bs 1.0 5e-4) in
  let de1 = abs_float ((e1 +. 0.25) /. 0.25) and 
      de2 = abs_float ((e2 +. 0.25) /. 0.25) in 
  let r = de1 /. de2 in 
    assert_bool (Printf.sprintf "energy scaling = %g (should be about 4)" r)
      ((r < (sqrt 32.0)) && (r > (sqrt 8.0)))

let tests = "leapfrog.ml tests" >:::
  ["const-timestep leapfrog energy scaling" >:: test_const_dt_energy;
   "leapfrog energy error scaling" >:: test_leap_energy]
