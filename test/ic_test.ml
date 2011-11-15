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

let test_standard_units () = 
  let bs = Array.init 1000
    (fun _ -> 
      B.make 0.0 1.0 
        (Array.init 3 (fun _ -> Random.float 1.0))
        (Array.init 3 (fun _ -> Random.float 0.01))) in
  let bs = Ic.adjust_frame bs in 
  let eini = E.energy bs in 
  let bs' = Ic.rescale_to_standard_units bs in 
  let efin' = E.energy bs' and 
      efin = E.energy bs and
      mtot = A.total_mass bs' in 
    assert_equal_float ~msg:"bs changed" efin eini;
    assert_equal_float ~msg:"bs' not in correct units" (-0.25) efin';
    assert_equal_float ~msg:"bs' total_mass wrong" 1.0 mtot

let test_shift_system () = 
  let bs = Ic.make_plummer 100 in 
  let v0 = Array.init 3 (fun _ -> Random.float 1.0 -. 0.5) and 
      r0 = Array.init 3 (fun _ -> Random.float 1.0 -. 0.5) in 
  let bs' = Ic.shift_system bs r0 v0 in 
  let com' = A.center_of_mass bs' and 
      ptot' = A.total_momentum bs' in 
    assert_equal_float_array com' r0;
    assert_equal_float_array ptot' v0

let test_generate_binary () = 
  for i = 0 to 100 do 
    let m1 = Random.float 1.0 and 
        m2 = Random.float 1.0 and 
        elts = Ic.random_elements 0.1 100.0 in 
    let bs = Ic.generate_binary m1 m2 elts in 
      assert_bool "Not in COM frame" ((Base.norm (A.center_of_mass bs)) < 1e-8);
      match bs with 
        | [|b1; b2|] -> 
          let e = A.binary_energy b1 b2 in 
          let a = elts.Ic.a in
          let ecc = A.eccentricity b1 b2 and 
              i = A.inclination b1 b2 in 
            assert_equal_float ~msg:"sma doesn't match energy" (~-.m1*.m2/.(2.0*.a)) e;
            assert_equal_float ~msg:"eccentricity doesn't match" ecc elts.Ic.e;
            assert_equal_float ~msg:"inclination doesn't match" i elts.Ic.i
        | _ -> raise (Failure "generate binary didn't produce two bodies")
  done

let tests = "ic.ml tests" >:::
  ["plummer model energy test" >:: test_plummer_energies;
   "make_{hot,cold}_spherical energy test" >:: test_spherical_energies;
   "plummer lagrange radii" >:: test_plummer_lagrange_radii;
   "cold_spherical lagrange radii" >:: test_cold_spherical_lagrange_radii;
   "hot_spherical lagrange radii" >:: test_hot_spherical_lagrange_radii;
   "standard units test" >:: test_standard_units;
   "shift_system test" >:: test_shift_system;
   "generate_binary test" >:: test_generate_binary]
