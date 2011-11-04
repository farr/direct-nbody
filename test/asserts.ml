open OUnit

let assert_equal_float ?(epsabs = 1e-8) ?(epsrel = 1e-8) ?(msg = "") x y = 
  assert_equal ~cmp:(cmp_float ~epsabs:epsabs ~epsrel:epsrel) ~printer:string_of_float ~msg:msg
    x y

let assert_equal_float_array ?(epsabs = 1e-8) ?(epsrel = 1e-8) ?(msg = "") xs ys = 
  let n = Array.length xs in 
    assert_equal ~printer:string_of_int ~msg:(msg ^ ": array lengths differ") n (Array.length ys);
    for i = 0 to n - 1 do 
      assert_equal_float ~epsabs:epsabs ~epsrel:epsrel ~msg:(msg ^ ": component " ^ (string_of_int i) ^ " differs") xs.(i) ys.(i)
    done
