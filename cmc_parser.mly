%token <float> FLOAT
%token EOL EOF
%start main
%type <float array array> main
%%
main:
float_arrays EOF { Array.of_list (List.rev $1) }
| float_arrays EOL EOF { Array.of_list (List.rev $1) }
;

float_arrays:
  float_array { [Array.of_list (List.rev $1)] }
| float_arrays EOL float_array { (Array.of_list (List.rev $3)) :: $1 }
;

float_array:
  FLOAT { [$1] }
| float_array FLOAT { $2 :: $1 }
;
