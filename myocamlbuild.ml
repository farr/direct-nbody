open Ocamlbuild_plugin

let oUnit_dir = "/Users/farr/Documents/code/ounit-1.0.3"

let _ = dispatch begin function 
  | After_rules -> 
      ocaml_lib ~extern:true ~dir:oUnit_dir "oUnit"
  | _ -> ()
end
