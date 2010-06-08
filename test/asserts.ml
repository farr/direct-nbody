open OUnit

let assert_equal_float ?(epsabs = 1e-8) ?(epsrel = 1e-8) ?(msg = "") x y = 
  assert_equal ~cmp:(cmp_float ~epsabs:epsabs ~epsrel:epsrel) ~printer:string_of_float ~msg:msg
    x y
