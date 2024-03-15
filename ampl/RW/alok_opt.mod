reset;
option solver ipopt;
# options ipopt_options 'nlp_scaling_method=gradient-based max_iter=1000 print_level=3';
set J :={ 1 ,2  };
set H :={ 1 ,2  };
set JH :={J,H};
var c >= 0 , <= 1.0;
var nh {JH}>=2 ;
param S2 {JH} ;
param Nd {JH} ;
param Y {J} ;

subject to ogrc {i in J}: 
  sum {k in H} (Nd[i,k]*Nd[i,k]*(1/nh[i,k]-1/Nd[i,k])*S2[i,k]) <= c*c*Y[i]*Y[i];
subject to ogrN0 {(i,k) in JH} : nh[i,k]>=2;
subject to ogrN {(i,k) in JH} : nh[i,k]<=Nd[i,k];
subject to  ogrnntot :  sum {(i,k) in JH} (nh[i,k])<=  1400  ;
minimize precyzja :      c ;
data;
param S2 := 
[ 1 , 1 ] 0.927338025672348
[ 1 , 2 ] 45.7014929370454
[ 2 , 1 ] 0.55133538255141
[ 2 , 2 ] 31.8778390978892
; 
param Nd := 
[ 1 , 1 ] 1000
[ 1 , 2 ] 1000
[ 2 , 1 ] 87
[ 2 , 2 ] 113
; 
param Y :=
1 16936.3869863486
2 1594.41022989484
; 
solve;
let {(i,k) in JH} nh[i,k]:=round(nh[i,k]);
printf {(i,k) in JH}: "%10d \n",nh[i,k] > "nh.txt"
;
quit;
