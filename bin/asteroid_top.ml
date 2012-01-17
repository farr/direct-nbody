(** Progam written to test what happens to the orbit of a pair
    of asteroids as their orbit about the sun undergoes Kozai
    oscillations. *)

open Base

module Ad = Leapfrog
module Ic = Ics.Make(Ad)
module An = Analysis.Make(Ad)
module E = An.E

let pi = 3.1415926535897932385

let mast = 1e-8
let mjup = 1e-3
let msun = 1.0

let aast = 1.0
let ajup = 10.0

let e0ast = 1e-3

let sepast = 1e-4

let i0jup = 80.0 *. pi /. 180.0
let i0ast = 0.0

let pjup = 2.0*.pi /. (sqrt ((msun +. mjup) /. ajup**3.0))
let past = 2.0*.pi /. (sqrt ((msun +. 2.0*.mast) /. aast**3.0))

let protast = 2.0 *. pi /. (sqrt (2.0*.mast /. sepast**3.0))

let pkozai = 2.0*.pjup**2.0/.(3.0*.pi*.past)*.(msun +. 2.0*.mast +. mjup) /. mjup

let dt = (min (min pjup past) protast) /. 200.0

let nout = int_of_float (floor (pkozai /. (1000.0*.dt)))
 
let step bs = 
  Array.iter (fun b -> Ad.drift (0.5*.dt) b) bs;
  for i = 0 to Array.length bs - 1 do
    let bi = bs.(i) in 
      for j = i+1 to Array.length bs - 1 do 
        let bj = bs.(j) in 
          Ad.kick infinity dt bi bj
      done
  done;
  Array.iter (fun b -> Ad.drift (0.5*.dt) b) bs

let afilter pred a = 
  Array.of_list (List.filter pred (Array.to_list a))

let angle l1 l2 = 
  acos ((dot l1 l2) /. (norm l1) /. (norm l2))

let mutual_inclination bs1 bs2 = 
  match bs1, bs2 with 
    | [|b11; b12|], [|b21; b22|] -> 
      let l1 = An.binary_angular_momentum b11 b12 and 
          l2 = An.binary_angular_momentum b21 b22 in 
        acos ((dot l1 l2) /. (norm l1) /. (norm l2))
    | _ -> 
      raise (Invalid_argument "mutual_inclination")

let asts = 
  let m = Random.float (2.0*.pi) in 
    Ic.generate_binary mast mast {Ic.m = m; a = sepast; e = e0ast; i = i0ast; capom = 0.0; omega = 0.0}

let jupsun = 
  Ic.generate_binary mjup msun {Ic.m = Random.float (2.0*.pi); a = ajup; i = i0jup; e = 0.0; capom = 0.0; omega = 0.0}

let (_,_,rast,vast) = 
  Ic.elements_to_rv msun (2.0*.mast) {Ic.m = Random.float (2.0*.pi); a = aast; i = i0ast; e = 0.0; capom = 0.0; omega = 0.0}

let bs = Ic.adjust_frame (Ic.combine_systems jupsun (Ic.shift_system asts rast vast))

let get_asts bs = afilter (fun b -> Ad.m b = mast) bs
let get_jupsun bs = afilter (fun b -> not (Ad.m b = mast)) bs
let get_sun bs = (afilter (fun b -> Ad.m b = msun) bs).(0)

let mutual_inclination b1 b2 b3 b4 = 
  let l1 = An.binary_angular_momentum b1 b2 and 
      l2 = An.binary_angular_momentum b3 b4 in 
    acos ((dot l1 l2) /. (norm l1) /. (norm l2))

let output bs = 
  let asts = get_asts bs and 
      jupsun = get_jupsun bs and 
      sun = get_sun bs in 
  let ast_summary = An.summary_body asts in 
  let aast = An.semi_major_axis asts.(0) asts.(1) and 
      east = An.eccentricity asts.(0) asts.(1) and 
      iast = An.inclination asts.(0) asts.(1) in 
  let asast = An.semi_major_axis ast_summary sun and 
      esast = An.eccentricity ast_summary sun and 
      isast = An.inclination ast_summary sun in 
  let ajup = An.semi_major_axis jupsun.(0) jupsun.(1) and 
      ejup = An.eccentricity jupsun.(0) jupsun.(1) and 
      ijup = An.inclination jupsun.(0) jupsun.(1) in 
  let imutsunasts = mutual_inclination sun ast_summary asts.(0) asts.(1) and
      imutjupast = mutual_inclination sun ast_summary jupsun.(0) jupsun.(1) in
    Printf.printf "%g %g %g %g %g %g %g %g %g %g %g %g\n%!"
      (Ad.t bs.(0)) aast east iast asast esast isast ajup ejup ijup imutsunasts imutjupast

let _ = 
  output bs;
  while true do 
    for i = 1 to nout do 
      step bs
    done;
    output bs
  done
