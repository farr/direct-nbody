open OUnit
open Asserts

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

let test_plummer_lagrange_radii () = 
  let bs = Ic.make_plummer 1000 in 
  let quarter = A.lagrange_radius bs 0.25 and 
      half = A.lagrange_radius bs 0.5 and 
      tquarter = A.lagrange_radius bs 0.75 in 
    assert_equal_float ~epsrel:0.1 0.4778 quarter;
    assert_equal_float ~epsrel:0.1 0.7686 half;
    assert_equal_float ~epsrel:0.1 1.2811 tquarter

let test_cold_spherical_lagrange_radii () = 
  let bs = Ic.make_cold_spherical 1000 and 
      r = 12.0 /. 5.0 in 
    assert_equal_float ~epsrel:0.1 (0.629961*.r) (A.lagrange_radius bs 0.25);
    assert_equal_float ~epsrel:0.1 (0.793701*.r) (A.lagrange_radius bs 0.5);
    assert_equal_float ~epsrel:0.1 (0.90856*.r) (A.lagrange_radius bs 0.75)

let test_hot_spherical_lagrange_radii () = 
  let bs = Ic.make_hot_spherical 1000 in 
    assert_equal_float ~epsrel:0.1 0.629961 (A.lagrange_radius bs 0.25);
    assert_equal_float ~epsrel:0.1 0.793701 (A.lagrange_radius bs 0.5);
    assert_equal_float ~epsrel:0.1 0.90856 (A.lagrange_radius bs 0.75)

let tests = "ic.ml tests" >:::
  ["plummer model energy test" >:: test_plummer_energies;
   "make_{hot,cold}_spherical energy test" >:: test_spherical_energies;
   "plummer lagrange radii" >:: test_plummer_lagrange_radii;
   "cold_spherical lagrange radii" >:: test_cold_spherical_lagrange_radii;
   "hot_spherical lagrange radii" >:: test_hot_spherical_lagrange_radii]
