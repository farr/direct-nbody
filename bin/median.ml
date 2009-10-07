(** Computes the median of a list of samples, plus the "one-sigma"
    interval about the median using the bootstrap method. *)

let nsamp = ref 100
let low = ref 0.158655
let high = ref 0.841345

let _ = 
  Arg.parse
    [("-nsamp", Arg.Set_int nsamp, "number of bootstrap samples for errorbars");
     ("-low", Arg.Set_float low, "low errorbar cumulative probability");
     ("-high", Arg.Set_float high, "high errorbar cumulative probability")]
    (fun _ -> ())
    "median [OPTIONS ...]"

let fcmp (x : float) y = 
  if x < y then 
    -1
  else if x = y then 
    0 
  else
    1

let median (xs : float array) = 
  let n = Array.length xs in 
  let xs = Array.copy xs in 
  Array.fast_sort fcmp xs;
  if n mod 2 = 0 then 
    0.5*.(xs.(n/2) +. xs.(n/2-1))
  else
    xs.(n/2)

let bootstrap (xs : float array) = 
  let n = Array.length xs in 
  let bxs = Array.make n 0.0 in 
  for i = 0 to n - 1 do 
    bxs.(i) <- xs.(Random.int n)
  done;
  bxs

let round x = int_of_float (x +. 0.5)

let isnan x = match classify_float x with 
| FP_nan -> true
| _ -> false

let bootstrap_median xs = 
  let xs = 
    Array.of_list (List.filter (fun x -> not (isnan x)) (Array.to_list xs)) in
  let meds = Array.make !nsamp 0.0 in 
  for i = 0 to !nsamp - 1 do 
    meds.(i) <- median (bootstrap xs)
  done;
  Array.fast_sort fcmp meds;
  let ilow = round (!low*.(float_of_int (!nsamp - 1))) and 
      ihigh = round (!high*.(float_of_int (!nsamp - 1))) in 
  (median xs, meds.(ilow), meds.(ihigh))

let _ = 
  let nums = 
    Array.of_list
      (let nums = ref [] in 
      try 
        while true do 
          nums := (Scanf.fscanf stdin "%s\n" float_of_string) :: !nums
        done;
        !nums
      with 
      | End_of_file -> 
          !nums) in 
  let (med, low, high) = bootstrap_median nums in 
  Printf.printf "%g %g %g\n" med low high
