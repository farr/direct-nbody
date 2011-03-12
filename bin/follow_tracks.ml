module A = Advancer.A
module Ic = Ics.Make(A)
module E = Energy.Make(A)

let nbs = ref 16
let dt = ref 0.1
let eta = ref 1e-5
let de = ref 1e-3

let options = 
  [("-n", Arg.Set_int nbs,
    "N number of bodies");
   ("-dt", Arg.Set_float dt,
    "DT timestep for output");
   ("-eta", Arg.Set_float eta,
    "ETA safety factor for timestep adjustment");
   ("-de", Arg.Set_float de,
    "DE maximum relative error in energy tolerated")]

let print_array v = 
  Array.iter (fun x -> Printf.printf "%20.15g " x) v

let _ = 
  Arg.parse 
    (Arg.align options)
    (fun _ -> ())
    "follow_tracks.native OPTIONS ...";
  Random.self_init ();
  A.follow_flag := true;
  let bs = Ic.rescale_to_standard_units (Ic.adjust_frame (Ic.make_plummer !nbs)) in 
  let filter bs = 
    Array.iter 
      (fun b -> 
        Printf.printf "%4d %20.15g %20.15g " b.A.id b.A.m b.A.t;
        print_array b.A.q;
        let v = Array.map (fun p -> p /. b.A.m) b.A.p in
          print_array v;
          Printf.printf "\n%!")
      bs;
    Printf.printf "\n";
    bs in
  let e0 = E.energy bs in
  let rec loop bs = 
    let new_bs = A.advance bs !dt !eta in
    let e = E.energy new_bs in 
    let delta_e = abs_float ((e -. e0) /. e0) in 
      if delta_e > !de then begin
        Printf.eprintf "Energy error limit exceeded!";
        exit 1
      end else begin
        loop new_bs
      end in
    loop (filter bs)
          
    
