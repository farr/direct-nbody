open OUnit

module B = Advancer.A
module Ic = Ics.Make(B)
module A = Analysis.Make(B)
module E = A.E

let test_plummer_energies () = 
  let bs = Ic.make_plummer 1000 in 
  let ke = E.total_kinetic_energy bs and 
      pe = E.total_potential_energy bs in 
  assert_equal ~cmp:(cmp_float ~epsrel:0.07) ke 0.25;
  assert_equal ~cmp:(cmp_float ~epsrel:0.07) pe (-0.5)

let tests = "ic.ml tests" >:::
  ["plummer model energy test" >:: test_plummer_energies]
