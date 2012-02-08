open OUnit

module B = Advancer.A
module Ic = Ics.Make(B)
module A = Analysis.Make(B)
module E = A.E
module Ad = Advancer.A

let test_energy_error () = 
  let n = 10 in 
  let bs = Ic.rescale_to_standard_units (Ic.make_plummer n) in 
  Base.set_eps (4.0/.(float_of_int n));
  let new_bs = Ad.advance bs 1.0 1e-4 in 
  let new_bs2 = Ad.advance bs 1.0 1e-3 in 
  let e0 = E.energy bs and 
      e1 = E.energy new_bs and 
      e2 = E.energy new_bs2 in 
  let de1 = abs_float (e1 -. e0) and 
      de2 = abs_float (e2 -. e0) in 
  let r = de2 /. de1 in 
  let r_exact = 10.0**(4.0/.3.0) in
  let ratio = r /. r_exact in 
    Base.set_eps 0.0;
  (* Ensure within 1/2 logarithmic order of magnitude *)
    assert_bool "ratio too small" (ratio > (sqrt 0.1)); 
    assert_bool "ratio too big" (ratio < (sqrt 10.0));
    assert_equal ~msg:"energy error too big" 
      ~cmp:(cmp_float ~epsabs:0.0 ~epsrel:1e-6)
      ~printer:string_of_float
      e1
      e0

(* Plummer model plus M = 1 central Plummer potential with eps = 0.1. *)
let test_external_potential () = 
  let n = 100 in 
  let bs = Ic.rescale_to_standard_units (Ic.make_plummer n) in 
    Base.set_eps (4.0 /. (float_of_int n));
  let v b = 
    let m = b.B.m and r = Base.norm (b.B.q) in 
    let r2 = r*.r +. 0.01 in (* r^2 + eps^2 with eps = 0.1 *)
    ~-.m /. (sqrt r2) in 
  let gv b = 
    let m = b.B.m and r = Base.norm b.B.q in 
    let r2 = r*.r +. 0.01 in 
    let r3 = r2*.(sqrt r2) in 
      Array.map (fun x -> m*.x/.r3) b.B.q in 
  let bs2 = Ad.advance ~extpot:gv bs 1.0 1e-2 in 
  let e0 = E.energy bs +. (E.total_external_energy v bs) and 
      e1 = E.energy bs2 +. (E.total_external_energy v bs2) in 
    Base.set_eps 0.0;
  let de = abs_float ((e1 -. e0) /. e0) in 
    assert_bool (Printf.sprintf "energy error too large: %g" de) (de < 1e-2)

let tests = "advancer.ml tests" >:::
  ["energy error test" >:: test_energy_error;
   "external potential test" >:: test_external_potential]
