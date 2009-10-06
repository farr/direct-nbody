open OUnit

module B = Advancer.A
module Ic = Ics.Make(B)
module A = Analysis.Make(B)
module E = A.E

let test_e maker n ke pe = 
  let bs = maker n in 
  let bke = E.total_kinetic_energy bs and 
      bpe = E.total_potential_energy bs in 
  let sigma = 1.0 /. (sqrt (float_of_int n)) in 
  let cmp = cmp_float ~epsrel:(3.0*.sigma) ~epsabs:(3.0*.sigma) in 
  assert_equal ~cmp:cmp ~printer:string_of_float bke ke;
  assert_equal ~cmp:cmp ~printer:string_of_float bpe pe

let test_plummer_energies () = 
  test_e Ic.make_plummer 1000 0.25 (-0.5)

let test_spherical_energies () = 
  test_e Ic.make_cold_spherical 1000 0.0 (-0.25);
  test_e Ic.make_hot_spherical 1000 0.25 (-0.5)

let tests = "ic.ml tests" >:::
  ["plummer model energy test" >:: test_plummer_energies;
   "make_{hot,cold}_spherical energy test" >:: test_spherical_energies]
