# Alokacja optymalna w domenach w schemacie warstwowym z uzyciem AMPL

AMPL_sub<-function(n=n,data=din,J="wo_su",H="h",S2h="S2h",Nh="Nh",Yhat="Yhat")
  # alokacja optymalna za pomoca programu AMPL
  # wariant z podpopulacjami
{
  
  nJ<-length(unique(data[[J]]))
  nH<-length(unique(data[[H]]))
  S2<-matrix(data[[S2h]],nJ,nH,byrow=TRUE)
  Nd<-matrix(data[[Nh]],nJ,nH,byrow=TRUE)
  Y<-matrix(data[[Yhat]],nJ,nH,byrow=TRUE)
  Y<-Y[,1]
  
  dir_old <- getwd()
  setwd('E:/amplcml')
  plik<-'alok_opt.mod'
  
  txt.mod<-function(a)
    write.table(a,plik,append=TRUE,quote=FALSE,row.names=FALSE,col.names=FALSE)
  
  if (file.exists(plik)==TRUE) file.remove(plik)
  
  txt.mod("reset;")
  txt.mod("option solver ipopt;")
  txt.mod("# options ipopt_options 'nlp_scaling_method=gradient-based max_iter=1000 print_level=3';")
  #txt.mod('option solver sdonlp2;')
  #txt.mod("options minos_options 'Major_iterations=50000 Superbasics_limit=1000 feasibility_tolerance=1.0e-5';")
  
  subs<-noquote(unique(data[[J]]))
  subs<-paste(subs,"",collapse=",")
  
  hs<-noquote(unique(data[[H]]))
  hs<-paste(hs,"",collapse=",")
  
  txt.mod(paste('set J :={',subs,'};'))
  txt.mod(paste('set H :={',hs,'};'))
  txt.mod('set JH :={J,H};')
  
  
  txt.mod('var c >= 0 , <= 1.0;');
  #txt.mod('var nh {JH}>=2, integer ;');
  txt.mod('var nh {JH}>=2 ;');
  txt.mod('param S2 {JH} ;')
  txt.mod('param Nd {JH} ;')
  txt.mod('param Y {J} ;')  
  
  txt.mod('')
  
  txt.mod('subject to ogrc {i in J}: ')  
  txt.mod('  sum {k in H} (Nd[i,k]*Nd[i,k]*(1/nh[i,k]-1/Nd[i,k])*S2[i,k]) <= c*c*Y[i]*Y[i];')
  
  txt.mod('subject to ogrN0 {(i,k) in JH} : nh[i,k]>=2;')
  txt.mod('subject to ogrN {(i,k) in JH} : nh[i,k]<=Nd[i,k];')
  
  txt.mod(paste('subject to  ogrnntot :  sum {(i,k) in JH} (nh[i,k])<= ',n,' ;'))
  
  txt.mod('minimize precyzja :      c ;');
  
  
  txt.mod('data;')
  txt.mod('param S2 := ')
  
  for (i in 1:nrow(data)) txt.mod(paste('[',data[i,J],',',data[i,H],']',data[i,S2h]))
  txt.mod('; ')
  
  
  txt.mod('param Nd := ')
  for (i in 1:nrow(data)) txt.mod(paste('[',data[i,J],',',data[i,H],']',data[i,Nh]))
  txt.mod('; ')
  
  
  txt.mod('param Y :=')
  ii=0;
  for (i in unique(data[[J]])) 
  {
    ii=ii+1
    txt.mod(paste(i,Y[ii]))  
  }
  txt.mod('; ')
  
  txt.mod('solve;')
  txt.mod('let {(i,k) in JH} nh[i,k]:=round(nh[i,k]);')
  #txt.mod(paste('display nh > ','"','nh.txt','"',sep=""))
  #txt.mod("let {(i,k) in JH} nh[i,k]:="); 
  #txt.mod("  (if Uniform01()>=ceil(nh[i,k])-nh[i,k] then ceil(nh[i,k]) ") # korekta <=;
  #txt.mod("  else floor(nh[i,k])); ")
  txt.mod( paste('printf {(i,k) in JH}: ','"','%10d \\n','"',',nh[i,k] > ','"nh.txt','"',sep="") )
  txt.mod(';')
  txt.mod('quit;')
  
  system('ampl.exe alok_opt.mod',invisible=TRUE)
  nh<-read.table('nh.txt',nrows=nrow(data))
  nh<-nh[,1]
  
  setwd(dir_old)
  return(nh)
}
