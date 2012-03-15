%token <float> FLOAT
%token EOL
%start main
%type <float array array> main
%%
main:
float_arrays { Array.of_list (List.rev $1) }
;

float_arrays:
  float_arrays float_array EOL { (Array.of_list (List.rev $2)) :: $1 }
;

float_array:
  float_array FLOAT { $2 :: $1 }
;
