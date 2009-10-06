(*  integrator.ml: Drives advancers in a long-time integration.
    Copyright (C) 2006--2008 Will M. Farr <farr@mit.edu>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*)

open Body;;
open Base;;
open Advancer;;

module type INTEGRATOR = 
  sig
    type b 
    exception Eg_err of b array
    val integrate : ?ci : float -> ?sf : float -> 
      ?max_err : float -> ?step_err : float -> 
	?filter : (b array -> b array) -> 
          ?max_dt : (b array -> float) -> 
            ?eg : float -> b array -> float 
	    -> b array
    val constant_sf_integrate : ?ci : float -> ?sf : float -> 
      ?filter : (b array -> b array) -> 
        ?max_dt : (b array -> float) -> 
          ?max_eg_err : float -> 
        b array -> float -> b array
  end

module Make(B : BODY)(A : (ADVANCER with type b = B.b)) : (INTEGRATOR with type b = A.b) = 
  struct 
    type b = A.b

    module Eg = Energy.Make(B)

    exception Eg_err of b array

    let adjust_sf sf err target = 
      let try_sf = sf *. (target/.err) in 
      if try_sf > sf *. 2.0 then 
	sf *. 2.0
      else if try_sf < sf /. 5.0 then 
	sf/.5.0
      else
	try_sf

    let integrate ?ci ?(sf = 0.1) 
	?(max_err = 1e-3) ?(step_err = 1e-6) 
	?(filter = fun (bs : b array) -> bs)
        ?(max_dt = fun (bs : b array) -> infinity)
	?eg bs dt = 
      let ci = match ci with 
      | None -> dt
      | Some(t) -> t in 
      let eg = match eg with 
      | None -> Eg.energy bs
      | Some(e) -> e in 
      let stop_time = (B.t bs.(0)) +. dt in 
      let rec loop cbs last_eg sf = 
	if B.t cbs.(0) >= stop_time then 
	  cbs
	else
          let dt_max = max_dt bs in 
	  let nbs = A.advance cbs (min ci dt_max) sf in 
	  let neg = Eg.energy nbs in 
	  let err = relative_error neg eg in 
	  let serr = relative_error neg last_eg in 
	  if serr > 5.0*.step_err then 
	    begin
	      let new_sf = adjust_sf sf serr step_err in 
	      Printf.fprintf stderr 
		"Retrying step at time %g due to step error %g (new sf: %g)\n"
		(B.t cbs.(0)) serr new_sf;
	      flush stderr;
	      loop cbs last_eg new_sf
	    end
	  else if err > max_err then 
	    raise (Eg_err cbs)
	  else
	    begin
	      let new_sf = adjust_sf sf serr step_err in 
	      Printf.fprintf stderr
		"Successful step to %g with step error %g (new sf: %g)\n"
		(B.t nbs.(0)) serr new_sf;
	      flush stderr;
	      loop (filter nbs) neg new_sf
	    end in 
      loop bs eg sf

    let constant_sf_integrate 
        ?(ci = 1.0) 
        ?(sf = 1e-4) 
        ?(filter = fun (bs : b array) -> bs)
        ?(max_dt = fun (bs : b array) -> infinity)
        ?(max_eg_err = 1e-5)
        bs dt = 
      let e = Eg.energy bs and 
          stop_t = dt +. (B.t bs.(0)) in 
      let rec loop bs = 
        let dt_max = max_dt bs in 
        let new_bs = A.advance bs (min dt_max ci) sf in 
        let new_e = Eg.energy new_bs and 
            new_t = B.t new_bs.(0) in 
        let eerr = abs_float ((new_e -. e)/.e) in 
        if eerr > max_eg_err then 
          raise (Eg_err bs)
        else begin
          Printf.fprintf stderr 
            "Just stepped to %g with energy error %g\n%!"
            new_t eerr;
          if new_t >= stop_t then 
            filter new_bs
          else
            loop (filter new_bs)
        end in 
      loop bs
        
  end
