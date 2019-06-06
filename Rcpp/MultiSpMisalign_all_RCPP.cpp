# include <RcppArmadillo.h>
# include <RcppArmadilloExtensions/sample.h>
# include <omp.h>

// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

double logitcpp(double theta,double a,double b){
  return log((theta-a)/(b-theta));
}

double logitinvcpp(double z,double a,double b){
  return b-(b-a)/(1+exp(z));
} 


double spCor(double D, double phi){
  return exp(-1.0*phi*D);
}

// [[Rcpp::export]]
arma::mat mkCovcpp2(mat A, colvec rho, mat D, ivec n, colvec psi){// input is distance matrix
  mat K=A*A.t();
  mat V=D;
  int q = A.n_rows;
  int i,j,jj,kk,k;
  
  for(i=0,j=0; i<q;++i){
    if(n[i]!=0) {
      V(span(j,j+n[i]-1),span(j,j+n[i]-1))=K(i,i)*exp(-rho[i]*V(span(j,j+n[i]-1),span(j,j+n[i]-1)));
      (V(span(j,j+n[i]-1),span(j,j+n[i]-1))).diag() += psi[i];
    }
    j += n[i];
  }
  for(i=0,jj=0;i<q-1;++i){
    
    kk=sum(n(span(0,i)));
    for(j=i+1;j<q; ++j){
      
      if(n[i]!=0 and n[j]!=0){
        arma::mat tempV = V(span(jj,jj+n[i]-1),span(kk,kk+n[j]-1));
        arma::mat temp = zeros(n[i],n[j]);
        for(k=0;k<q;++k){
          temp += A(i,k)*A(j,k)*exp(-rho[k]*tempV);
        }
        V(span(jj,jj+n[i]-1),span(kk,kk+n[j]-1)) = temp;
        V(span(kk,kk+n[j]-1),span(jj,jj+n[i]-1)) = temp.t();
      }
      // Rcout << "j:"<<kk<<"k:"<<kk+n[j]-1 <<"\n";
      kk += n[j];
    }
    // Rcout << "k:"<<jj<<"k:"<<jj+n[i]-1 <<"\n";
    jj += n[i];
  }
  return V;
}

// [[Rcpp::export]]
arma::mat mkCovcpp2_par(mat A, colvec rho, mat C, mat coordsD, ivec n, colvec psi,int N){// input is distance matrix
  
  omp_set_num_threads(2);
  int sl, sk, k,l,kk,jj,ii;
  sl = sk = 0;
  int m = 6;
  
  for(k = 0; k < m; k++){
    sl = 0;
    for(l = 0; l < m; l++){
      
#pragma omp parallel for private(kk,jj,ii)
      for(kk = 0; kk < n[k]; kk++){
        for(jj = 0; jj < n[l]; jj++){
          C[(sl+jj)*N+(sk+kk)] = 0.0;
          for(ii = 0; ii < m; ii++){
            // printf("%d ",k+m*ii);
            C[(sl+jj)*N+(sk+kk)] += A[k+m*ii]*A[l+m*ii]*exp(-1.0*coordsD[(sl+jj)*N+(sk+kk)]*rho[ii]);
          }
        }
      }
      sl += n[l];
    }
    sk += n[k];
  }
  
    sl = 0;
    for(l = 0; l < m; l++){
      for(k = 0; k < n[l]; k++){
        C[(sl+k)*N+(sl+k)] += psi[l];
      }
      sl += n[l];
    }
    return C;
}

// [[Rcpp::export]]
List MultiSpMisaligncpp_all(List Ylist_r, List Xlist_r, List Dlist_r, NumericVector Nlist_r,IntegerMatrix nlist,
                            int nt, int q, int p, NumericVector Aidx_r,NumericVector mAidx_r,NumericVector AMidx_r,
                            NumericVector A_step, double rho_step, double psi_step,
                            colvec rho_init, colvec psi_init, colvec beta_init, colvec A_init,
                            double sd_beta, double a_psi, double b_psi, mat SK,
                            int IW_df,double a_rho,double b_rho,
                            int iters, int n_report,bool fix_psi, int report){ // 2/15/2018 changed lognormal to normal prior for diag(A)
  ivec Nlist = as<ivec>(Nlist_r);
  uvec Aidx  = as<uvec>(Aidx_r);
  uvec mAidx = as<uvec>(mAidx_r);
  uvec AMidx = as<uvec>(AMidx_r);
  //set up:
  int nltr = q*(q+1)/2;
  int t;
  int pq = beta_init.n_elem;
  //initial values:
  field<mat> Xlist(nt,1);
  field<vec> Ylist(nt,1);
  field<mat> Dlist(nt,1);
  for(t = 0; t < nt; ++t){
    Xlist(t,0) = as<mat>(Xlist_r[t]);
    Ylist(t,0) = as<vec>(Ylist_r[t]);
    Dlist(t,0) = as<mat>(Dlist_r[t]);
  }
  colvec beta = beta_init;
  // double betafoo;
  colvec psi  = psi_init;
  colvec psi_can(q),rhocan(q);
  colvec rho  = rho_init;
  //make covariance matrix
  // form A, qxq lower triangle matrix
  colvec Aeach(nltr),Aeach_can(nltr),Aeach_samp(nltr);
  Aeach = A_init;
  // Aeach.elem(Aidx) = ones(q)*1.0;
  mat A=zeros(q,q);
  mat A_can=zeros(q,q);
  A.elem(AMidx) = Aeach;
  colvec Anorm = zeros(Aeach.size());
  // #save samples:
  mat keep_beta(pq,iters), keep_A(nltr,iters), keep_psi(q,iters);
  colvec accept=zeros(nltr);
  colvec accept_psi=zeros(q);
  colvec keep_rho=zeros(iters);
  colvec accept_rho=zeros(1);
  //for beta;
  mat XtV_spinvX(pq,pq); 
  colvec XtV_spinvY(pq);
  field<colvec> rlist(nt,1);
  //for psi
  double rtV_spinv_r,logDetV_sp,logpost,logpost_can,logMHRatio;
  field<mat> V_sp_chol(t,1);
  colvec temp(1);
  // for A
  colvec Atemp(nltr);
  // for rho
  double rho_t;
  // #Start MCMC:
  for(int ii = 0; ii < iters; ++ii){
    // Rcout <<"iter"<<ii<<"\n";
    //   ####################################################
    //   ###### update regression coefs #####
    //   ####################################################
    
    XtV_spinvX=zeros(pq,pq);
    XtV_spinvY=zeros(pq);
    
    for(t = 0; t < nt; ++t){
      mat V_sp=zeros(Nlist[t],Nlist[t]);
      V_sp = mkCovcpp2(A,rho,Dlist_r[t],nlist.row(t),psi);
      V_sp = chol(V_sp,"lower");
      mat V_spinvX = solve(trimatu(V_sp.t()),solve(trimatl(V_sp),Xlist(t,0)));
      colvec V_spinvY = solve(trimatu(V_sp.t()),solve(trimatl(V_sp),Ylist(t,0)));
      
      XtV_spinvX += trans(Xlist(t,0))*V_spinvX;
      XtV_spinvY += trans(Xlist(t,0))*V_spinvY;
      V_sp_chol(t,0) = V_sp;
    }
    XtV_spinvX.diag() += 1/pow(sd_beta,2);
    XtV_spinvX = 0.5*XtV_spinvX + 0.5*XtV_spinvX.t(); //added 6/4 doesn't change the results, but getride of warning
    mat cholv = chol(XtV_spinvX,"lower");
    colvec M = solve(trimatu(cholv.t()),solve(trimatl(cholv),XtV_spinvY));
    // colvec M = solve(XtV_spinvX,XtV_spinvY);
    // NumericVector ranz = rnorm(p*q);
    // beta = M+solve(trimatu(cholv.t()),as<colvec>(ranz));
    // betafoo = sum(M+solve(trimatu(cholv.t()),as<colvec>(ranz)));
    
    beta = M+solve(trimatu(cholv.t()),as<colvec>(rnorm(pq)));
    // betafoo = sum(M+solve(trimatu(cholv.t()),as<colvec>(rnorm(p*q))));
    // Rcout<<"ranz"<<ranz<< "\n";
    // Rcout<<"sizebeta"<<sum(beta)<< "\n";
    // Rcout<<"sizebeta"<<betafoo<< "\n";
    
    for(t = 0; t < nt; ++t)  rlist(t,0) = Ylist(t,0)-Xlist(t,0)*beta;
    
    //   ####################################################
    //   ###### update the variances #####
    //   ####################################################:
    
    rtV_spinv_r = logDetV_sp = 0.0;
    for(t = 0; t < nt; ++t){
      // # V_sp = mkCov(A,rho,coords,psi); V_sp_chol = t(chol(V_sp)) # done above
      temp = (rlist(t)).t()*solve(trimatu((V_sp_chol(t,0)).t()),solve(trimatl(V_sp_chol(t,0)),rlist(t,0)));
      rtV_spinv_r += temp[0];
      logDetV_sp += sum(2.0*log((V_sp_chol(t,0)).diag()));
    }
    // printf("%f ",logDetV_sp);
    logpost = sum(-(1+a_psi)*log(psi)-b_psi/psi)+sum(log(psi))-0.5*logDetV_sp-0.5*rtV_spinv_r;
    
    psi_can = exp(psi_step*as<vec>(rnorm(q))+log(psi));
    rtV_spinv_r = logDetV_sp = 0.0;
    for(t = 0; t < nt; ++t){
      mat V_sp =  mkCovcpp2(A,rho,Dlist_r[t],nlist.row(t),psi_can);
      V_sp = chol(V_sp,"lower");
      
      temp = (rlist(t)).t()*solve(trimatu(V_sp.t()),solve(trimatl(V_sp),rlist(t,0)));
      rtV_spinv_r += temp[0];
      logDetV_sp += sum(2.0*log(V_sp.diag()));
    }
    // printf("%f ",rtV_spinv_r);
    // printf("%f ",logDetV_sp);
    
    logpost_can = sum(-(1+a_psi)*log(psi_can)-b_psi/psi_can)+sum(log(psi_can))-0.5*logDetV_sp-0.5*rtV_spinv_r;
    
    logMHRatio = logpost_can - logpost;
    if( log(runif(1)[0]) <= logMHRatio){
      logpost = logpost_can;
      accept_psi[0] = accept_psi[0] +1;
    }else{
      psi_can=psi;
    }
    
    psi = psi_can;
    
    //       # A
    Anorm = A_step*rnorm(Aeach.size()); 
    // Rcout<<"Anorm"<< Anorm <<"\n";
    
    Aeach_samp.elem(Aidx) = exp(Anorm.elem(Aidx)+log(Aeach.elem(Aidx)));
    Aeach_samp.elem(mAidx) = Anorm.elem(mAidx) + Aeach.elem(mAidx);
    
    A.elem(AMidx) = Aeach;
    // logDetK = sum(2.0*log(diag(A))); SKtrace = sum(diag(backsolve(t(A),forwardsolve(A,SK))))
    rtV_spinv_r =  logDetV_sp = 0;
    for(t = 0; t < nt; ++t){
      mat V_sp =  mkCovcpp2(A,rho,Dlist_r[t],nlist.row(t),psi);
      V_sp = chol(V_sp,"lower");
      temp = (rlist(t)).t()*solve(trimatu(V_sp.t()),solve(trimatl(V_sp),rlist(t,0)));
      rtV_spinv_r += temp[0];
      logDetV_sp += sum(2.0*log(V_sp.diag()));
    }
    
    Atemp.elem(mAidx) =  Aeach.elem(mAidx);
    Atemp.elem(Aidx)  =  Aeach.elem(Aidx);
    // Rcout<<"Atemp"<< sum(log(A.diag())) -0.5*logDetV_sp-0.5*rtV_spinv_r-1.0/100.0*(Atemp.t()*Atemp)<<"\n";
    // colvec foo = NumericVector::create(Atemp.elem(Aidx));
    
    logpost = sum(log(A.diag()))-0.5*logDetV_sp-0.5*rtV_spinv_r-1.0/200.0*as_scalar(Atemp.t()*Atemp);
    Aeach_can = Aeach;
    Aeach_can = Aeach_samp;
    A_can.elem(AMidx) = Aeach_can;
    
    rtV_spinv_r =  logDetV_sp = 0;
    for(t = 0; t < nt; ++t){
      mat V_sp =  mkCovcpp2(A_can,rho,Dlist_r[t],nlist.row(t),psi);
      V_sp = chol(V_sp,"lower");
      temp = (rlist(t)).t()*solve(trimatu(V_sp.t()),solve(trimatl(V_sp),rlist(t,0)));
      rtV_spinv_r += temp[0];
      logDetV_sp += sum(2.0*log(V_sp.diag()));
    }
    
    Atemp.elem(mAidx) =  Aeach_can.elem(mAidx);
    Atemp.elem(Aidx)  =  Aeach_can.elem(Aidx);
    
    logpost_can = sum(log(A_can.diag()))-0.5*logDetV_sp-0.5*rtV_spinv_r-1.0/200.0*as_scalar(Atemp.t()*Atemp);
    logMHRatio = logpost_can - logpost;
    
    if(log(runif(1)[0]) <= logMHRatio){
      accept = accept+1;
      logpost = logpost_can;
    }else{
      Aeach_can = Aeach;
    }
    
    Aeach = Aeach_can;
    A.elem(AMidx) = Aeach;
    
    // # update rho
    rho_t = logitcpp(rho[0],a_rho,b_rho);
    rhocan = rep(logitinvcpp(rnorm(1,rho_t,0.1)[0],a_rho,b_rho),q);
    
    // # compute logpost
    rtV_spinv_r =  logDetV_sp = 0;
    for(t = 0; t < nt; ++t){
      
      mat V_sp =  mkCovcpp2(A,rho,Dlist_r[t],nlist.row(t),psi);
      V_sp = chol(V_sp,"lower");
      temp = (rlist(t)).t()*solve(trimatu(V_sp.t()),solve(trimatl(V_sp),rlist(t,0)));
      rtV_spinv_r += temp[0];
      logDetV_sp += sum(2.0*log(V_sp.diag()));
      
    }
    
    logpost = -0.5*logDetV_sp-0.5*rtV_spinv_r + log(rho[0] - a_rho) + log(b_rho - rho[0]);
    
    // # compute logpost_can
    rtV_spinv_r =  logDetV_sp = 0;
    for(t = 0; t < nt; ++t){
      
      mat V_sp =  mkCovcpp2(A,rhocan,Dlist_r[t],nlist.row(t),psi);
      V_sp = chol(V_sp,"lower");
      
      temp = (rlist(t)).t()*solve(trimatu(V_sp.t()),solve(trimatl(V_sp),rlist(t,0)));
      rtV_spinv_r += temp[0];
      logDetV_sp += sum(2.0*log(V_sp.diag()));
    }
    logpost_can = -0.5*logDetV_sp-0.5*rtV_spinv_r + log(rhocan[0] - a_rho) + log(b_rho - rhocan[0]);
    
    logMHRatio = logpost_can - logpost;
    
    if(log(runif(1)[0]) <= logMHRatio){
      accept_rho = accept_rho+1;
      rho = rhocan;
    }
    
    //       #save samples
    keep_beta.col(ii) = beta;
    keep_A.col(ii)    = Aeach;
    keep_psi.col(ii)  = psi;
    keep_rho[ii]  = rho[0];
    // //   //       # dev[i]<--2.0*sum(dnorm(r,0,1/sqrt(taue),log=T))
    // //   //       # if(spatial & i>burn){
    // //   //       #   theta.mn<-theta.mn+theta/(runs-burn)
    // //   //       #   theta.var<-theta.var+theta*theta/(runs-burn)
    // //   //       # }
    // //   //       #Display current iteration:
    if(ii%report==0){
      Rcout<<"--------------"<<ii<<"iterations-------------" << "\n";
      Rcout<<"acc:"<< accept/ii << accept_rho/ii << "\n";
    }
  }
  
  //     # theta.var<-theta.var-theta.mn^2
  //     # beta.mn<-apply(keep.beta[burn:runs,],2,mean)
  //     # sigma<-mean(1/sqrt(keep.taue[burn:runs]))
  //     # r<-Y-x%*%beta.mn
  //     # if(spatial){r<-r-P%*%theta.mn}
  //     # dhat<-sum(-2.0*dnorm(r,0,sigma,log=T))
  //     # dbar<-mean(dev[burn:runs])
  //     # pD<-dbar-dhat
  //     # DIC<-dbar+pD
  //     time = proc.time()-tick
  // return(list(beta=keep.beta, A = keep.A, psi = keep.psi, rho = keep.rho,acc.A=accept,acc.rho=accept_rho,time=time))
  // return keep_psi;
  return List::create(
    Named("beta") =keep_beta,
    Named("A") = keep_A,
    Named("psi") =keep_psi,
    Named("rho") = keep_rho,
    Named("accept.A") =accept,
    Named("accept.rho") = accept_rho);
}

// [[Rcpp::export]]
// parallel did not work for Rcpp
List MultiSpMisaligncpp_all_par(List Ylist_r, List Xlist_r, List Dlist_r, NumericVector Nlist_r,IntegerMatrix nlist,
                                int nt, int q, int p, NumericVector Aidx_r,NumericVector mAidx_r,NumericVector AMidx_r,
                                NumericVector A_step, double rho_step, double psi_step,
                                colvec rho_init, colvec psi_init, colvec beta_init, colvec A_init,
                                double sd_beta, double a_psi, double b_psi, mat SK,
                                int IW_df,double a_rho,double b_rho,
                                int iters, int n_report,bool fix_psi,int cores = 2){ // 2/15/2018 changed lognormal to normal prior for diag(A)
  omp_set_num_threads(cores);
  ivec Nlist = as<ivec>(Nlist_r);
  uvec Aidx  = as<uvec>(Aidx_r);
  uvec mAidx = as<uvec>(mAidx_r);
  uvec AMidx = as<uvec>(AMidx_r);
  //set up:
  int nltr = q*(q+1)/2;
  int t;
  //initial values:
  field<mat> Xlist(nt,1);
  field<vec> Ylist(nt,1);
  field<mat> Dlist(nt,1);
  for(t = 0; t < nt; ++t){
    Xlist(t,0) = as<mat>(Xlist_r[t]);
    Ylist(t,0) = as<vec>(Ylist_r[t]);
    Dlist(t,0) = as<mat>(Dlist_r[t]);
  }
  colvec beta = beta_init;

  colvec psi  = psi_init;
  colvec psi_can(q),rhocan(q);
  colvec rho  = rho_init;
  //make covariance matrix
  // form A, qxq lower triangle matrix
  colvec Aeach(nltr),Aeach_can(nltr),Aeach_samp(nltr);
  Aeach = A_init;
  // Aeach.elem(Aidx) = ones(q)*1.0;
  mat A=zeros(q,q);
  mat A_can=zeros(q,q);
  A.elem(AMidx) = Aeach;
  colvec Anorm = zeros(Aeach.size());
  // #save samples:
  mat keep_beta(p*q,iters), keep_A(nltr,iters), keep_psi(q,iters);
  colvec accept=zeros(nltr);
  colvec accept_psi=zeros(q);
  colvec keep_rho=zeros(iters);
  colvec accept_rho=zeros(1);
  //for beta;
  mat XtV_spinvX(p*q,p*q);
  colvec XtV_spinvY(p*q);
  field<colvec> rlist(nt,1);
  //for psi

  double rtV_spinv_r,logDetV_sp,logpost,logpost_can,logMHRatio;
  field<mat> V_sp_chol(t,1);
  colvec temp(1);
  // for A
  colvec Atemp(nltr);
  // for rho
  double rho_t;
  // #Start MCMC:
  for(int ii = 0; ii < iters; ++ii){
    // Rcout <<"iter"<<ii<<"\n";
    //   ####################################################
    //   ###### update regression coefs #####
    //   ####################################################

    XtV_spinvX=zeros(p*q,p*q);
    XtV_spinvY=zeros(p*q);
    
// #pragma omp parallel for private(t)
    for(t = 0; t < nt; ++t){
      mat V_sp=zeros(Nlist[t],Nlist[t]);
      V_sp = mkCovcpp2(A,rho,Dlist_r[t],nlist.row(t),psi);
      V_sp = chol(V_sp,"lower");
      mat V_spinvX = solve(trimatu(V_sp.t()),solve(trimatl(V_sp),Xlist(t,0)));
      colvec V_spinvY = solve(trimatu(V_sp.t()),solve(trimatl(V_sp),Ylist(t,0)));
      
// #pragma omp critical
      // {
        XtV_spinvX += trans(Xlist(t,0))*solve(trimatu(V_sp.t()),solve(trimatl(V_sp),Xlist(t,0)));
        XtV_spinvY += trans(Xlist(t,0))*solve(trimatu(V_sp.t()),solve(trimatl(V_sp),Ylist(t,0)));
      // }
      V_sp_chol(t,0) = V_sp;
    }

    XtV_spinvX.diag() += 1/pow(sd_beta,2);
    mat cholv = chol(XtV_spinvX,"lower");
    colvec M = solve(trimatu(cholv.t()),solve(trimatl(cholv),XtV_spinvY));
    // colvec M = solve(XtV_spinvX,XtV_spinvY);
    // NumericVector ranz = rnorm(p*q);
    beta = M+solve(trimatu(cholv.t()),as<colvec>(rnorm(p*q)));
    // Rcout<<"ranz"<<ranz<< "\n";
    //Rcout<<"sizebeta"<<size(beta)<< "\n";

    for(t = 0; t < nt; ++t)  rlist(t,0) = Ylist(t,0)-Xlist(t,0)*beta;

    //   ####################################################
    //   ###### update the variances #####
    //   ####################################################:

    rtV_spinv_r = logDetV_sp = 0.0;
#pragma omp parallel for private(temp,t)
    for(t = 0; t < nt; ++t){
      // # V_sp = mkCov(A,rho,coords,psi); V_sp_chol = t(chol(V_sp)) # done above
      temp = (rlist(t)).t()*solve(trimatu((V_sp_chol(t,0)).t()),solve(trimatl(V_sp_chol(t,0)),rlist(t,0)));
      
// #pragma omp critical
// {
      rtV_spinv_r += temp[0];
      logDetV_sp += sum(2.0*log((V_sp_chol(t,0)).diag()));
// }
    }
    // printf("%f ",logDetV_sp);
    logpost = sum(-(1+a_psi)*log(psi)-b_psi/psi)+sum(log(psi))-0.5*logDetV_sp-0.5*rtV_spinv_r;

    psi_can = exp(psi_step*as<vec>(rnorm(q))+log(psi));
    rtV_spinv_r = logDetV_sp = 0.0;

// #pragma omp parallel for private(temp,t)
    for(t = 0; t < nt; ++t){
      mat V_sp = mkCovcpp2(A,rho,Dlist_r[t],nlist.row(t),psi_can);
      V_sp = chol(V_sp,"lower");

      temp = (rlist(t)).t()*solve(trimatu(V_sp.t()),solve(trimatl(V_sp),rlist(t,0)));
      
// #pragma omp critical
// {
      rtV_spinv_r += temp[0];
      logDetV_sp += sum(2.0*log(V_sp.diag()));
// }
    }
    // printf("%f ",rtV_spinv_r);
    // printf("%f ",logDetV_sp);
    logpost_can = sum(-(1+a_psi)*log(psi_can)-b_psi/psi_can)+sum(log(psi_can))-0.5*logDetV_sp-0.5*rtV_spinv_r;

    logMHRatio = logpost_can - logpost;
    if( log(runif(1)[0]) <= logMHRatio){
      logpost = logpost_can;
      accept_psi[0] = accept_psi[0] +1;
    }else{
      psi_can=psi;
    }

    psi = psi_can;

    //       # A
    Anorm = A_step*rnorm(Aeach.size());
    Aeach_samp.elem(Aidx) = exp(Anorm.elem(Aidx)+log(Aeach.elem(Aidx)));
    Aeach_samp.elem(mAidx) = Anorm.elem(mAidx) + Aeach.elem(mAidx);

    A.elem(AMidx) = Aeach;
    // logDetK = sum(2.0*log(diag(A))); SKtrace = sum(diag(backsolve(t(A),forwardsolve(A,SK))))
    rtV_spinv_r =  logDetV_sp = 0;

// #pragma omp parallel for private(temp,t)
    for(t = 0; t < nt; ++t){

      mat V_sp =  mkCovcpp2(A,rho,Dlist_r[t],nlist.row(t),psi);
      V_sp = chol(V_sp,"lower");

      temp = (rlist(t)).t()*solve(trimatu(V_sp.t()),solve(trimatl(V_sp),rlist(t,0)));
      rtV_spinv_r += temp[0];
      logDetV_sp += sum(2.0*log(V_sp.diag()));
    }

    Atemp.elem(mAidx) =  Aeach.elem(mAidx);
    Atemp.elem(Aidx)  =  Aeach.elem(Aidx);



    logpost = sum(log(A.diag()))-0.5*logDetV_sp-0.5*rtV_spinv_r-1.0/200.0*as_scalar(Atemp.t()*Atemp);

    Aeach_can = Aeach;
    Aeach_can = Aeach_samp;
    A_can.elem(AMidx) = Aeach_can;

    rtV_spinv_r =  logDetV_sp = 0;
// #pragma omp parallel for private(temp,t)
    for(t = 0; t < nt; ++t){

      mat V_sp =  mkCovcpp2(A_can,rho,Dlist_r[t],nlist.row(t),psi);
      V_sp = chol(V_sp,"lower");

      temp = (rlist(t)).t()*solve(trimatu(V_sp.t()),solve(trimatl(V_sp),rlist(t,0)));
      rtV_spinv_r += temp[0];
      logDetV_sp += sum(2.0*log(V_sp.diag()));
    }

    Atemp.elem(mAidx) =  Aeach_can.elem(mAidx);
    Atemp.elem(Aidx)  =  Aeach_can.elem(Aidx);

    logpost_can = sum(log(A_can.diag()))-0.5*logDetV_sp-0.5*rtV_spinv_r-1.0/200.0*as_scalar(Atemp.t()*Atemp);

    logMHRatio = logpost_can - logpost;

    if(log(runif(1)[0]) <= logMHRatio){
      accept = accept+1;
      logpost = logpost_can;
    }else{
      Aeach_can = Aeach;
    }

    Aeach = Aeach_can;
    A.elem(AMidx) = Aeach;

    // # update rho
    rho_t = logitcpp(rho[0],a_rho,b_rho);
    rhocan = rep(logitinvcpp(rnorm(1,rho_t,0.1)[0],a_rho,b_rho),q);

    // # compute logpost
    rtV_spinv_r =  logDetV_sp = 0;
// #pragma omp parallel for private(temp,t)
    for(t = 0; t < nt; ++t){

      mat V_sp =  mkCovcpp2(A,rho,Dlist_r[t],nlist.row(t),psi);
      V_sp = chol(V_sp,"lower");

      temp = (rlist(t)).t()*solve(trimatu(V_sp.t()),solve(trimatl(V_sp),rlist(t,0)));
      rtV_spinv_r += temp[0];
      logDetV_sp += sum(2.0*log(V_sp.diag()));

    }

    logpost = -0.5*logDetV_sp-0.5*rtV_spinv_r + log(rho[0] - a_rho) + log(b_rho - rho[0]);

    // # compute logpost_can
    rtV_spinv_r =  logDetV_sp = 0;
// #pragma omp parallel for private(temp,t)
    for(t = 0; t < nt; ++t){
      mat V_sp =  mkCovcpp2(A,rhocan,Dlist_r[t],nlist.row(t),psi);
      V_sp = chol(V_sp,"lower");
      temp = (rlist(t)).t()*solve(trimatu(V_sp.t()),solve(trimatl(V_sp),rlist(t,0)));
      rtV_spinv_r += temp[0];
      logDetV_sp += sum(2.0*log(V_sp.diag()));
    }
    logpost_can = -0.5*logDetV_sp-0.5*rtV_spinv_r + log(rhocan[0] - a_rho) + log(b_rho - rhocan[0]);

    logMHRatio = logpost_can - logpost;

    if(log(runif(1)[0]) <= logMHRatio){
      accept_rho = accept_rho+1;
      rho = rhocan;
    }

    //       #save samples
    keep_beta.col(ii) = beta;
    keep_A.col(ii)    = Aeach;
    keep_psi.col(ii)  = psi;
    keep_rho[ii]  = rho[0];
    if(ii%1000==0){
    Rcout<<"--------------"<<ii<<"iterations-------------" << "\n";
    Rcout<<"acc:"<< accept/ii << accept_rho/ii << "\n";
    }
  }
  return List::create(
    Named("beta") =keep_beta,
    Named("A") = keep_A,
    Named("psi") =keep_psi,
    Named("rho") = keep_rho,
    Named("accept.A") =accept,
    Named("accept.rho") = accept_rho);
}

///*** R
// 
// # require(rbenchmark)
// # tick = proc.time()
// # set.seed(1)
// # out2<-MultiSpMisaligncpp_all_par(Ylist, Xlist, Dlist, Nlist,nlist,
// #                                  nt, q, p,(Aidx-1),(mAidx-1),(AMindx-1),
// #                                  A_step, rho_step, psi_step,
// #                                  rho_init,  psi_init,  beta_init, A_init,
// #                                  sd_beta,  a_psi,  b_psi,SK,
// #                                  IW_df, a_rho, b_rho,
// #                                  10,  n_report, fix_psi,2)
// # tock2 = proc.time()-tick
// #
// # tick = proc.time()
// # set.seed(1)
// # out1<-MultiSpMisaligncpp_all(Ylist, Xlist, Dlist, Nlist,nlist,
// #                              nt, q, p,(Aidx-1),(mAidx-1),(AMindx-1),
// #                              A_step, rho_step, psi_step,
// #                              rho_init,  psi_init,  beta_init, A_init,
// #                              sd_beta,  a_psi,  b_psi,SK,
// #                              IW_df, a_rho, b_rho,
// #                              10,  n_report, fix_psi)
// #
// # tock1 = proc.time()-tick
// 
// #
// # benchmark(
// # out1<-MultiSpMisaligncpp_all(Ylist, Xlist, Dlist, Nlist,nlist,
// #                              nt, q, p,(Aidx-1),(mAidx-1),(AMindx-1),
// #                              A_step, rho_step, psi_step,
// #                              rho_init,  psi_init,  beta_init, A_init,
// #                              sd_beta,  a_psi,  b_psi,SK,
// #                              IW_df, a_rho, b_rho,
// #                              3,  n_report, fix_psi),
// # out2<-MultiSpMisaligncpp_all_par(Ylist, Xlist, Dlist, Nlist,nlist,
// #                              nt, q, p,(Aidx-1),(mAidx-1),(AMindx-1),
// #                              A_step, rho_step, psi_step,
// #                              rho_init,  psi_init,  beta_init, A_init,
// #                              sd_beta,  a_psi,  b_psi,SK,
// #                              IW_df, a_rho, b_rho,
// #                              3,  n_report, fix_psi,2),
// # order="relative", replications=100)[,1:4]
// # //
// # // set.seed(1)
// # // tick = proc.time()
// # // out1=MultiSpMisalign_all(
// # //   y,  x,  coords, nt, q,p,
// # //   A.step = A_step,rho.step = rho_step,psi.step =psi_step,
// # //   rho.init =rho_init,psi.init= psi_init,fix_psi=F,
// # //   sd.beta=sd_beta,a.psi=a_psi,b.psi=b_psi,SK, IW.df = IW_df, a.rho = a_rho, b.rho =b_rho,iters = 10,n.report = 10)
// # // tock1 = proc.time()-tick
// # //
// */
