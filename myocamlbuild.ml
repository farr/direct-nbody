open Ocamlbuild_plugin

let oUnit_dir = "/Users/farr/Documents/code/oUnit"

let _ = dispatch begin function 
  | After_rules -> 
      ocaml_lib ~extern:true ~dir:oUnit_dir "oUnit";
      flag ["native"; "compile"]
        (S[A"-inline"; A"1000"; A"-unsafe"]);
      flag ["native"; "link"]
        (S[A"-inline"; A"1000"])
  | _ -> ()
end
