open OUnit

module B = Advancer.A
module Ic = Ics.Make(B)
module A = Analysis.Make(B)
module E = A.E
module Ad = Advancer.A

let test_energy_error () = 
  let bs = Ic.make_plummer 250 in 
  let new_bs = Ad.advance bs 1.0 1e-8 in 
  let new_bs2 = Ad.advance bs 1.0 1e-7 in 
  let e0 = E.energy bs and 
      e1 = E.energy new_bs and 
      e2 = E.energy new_bs2 in 
  let de1 = abs_float (e1 -. e0) and 
      de2 = abs_float (e2 -. e0) in 
  let r = de2 /. de1 in 
  assert_bool "ratio too small" (r > (sqrt 10.0));
  assert_bool "ratio too big" (r < (sqrt 1000.0));
  assert_equal ~msg:"energy error too big" 
    ~cmp:(cmp_float ~epsabs:0.0 ~epsrel:1e-6)
    ~printer:string_of_float
    e1
    e0

let tests = "advancer.ml tests" >:::
  ["energy error test" >:: test_energy_error]
