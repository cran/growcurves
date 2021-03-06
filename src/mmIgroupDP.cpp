#include "mmC.h"

using namespace Rcpp;
using namespace arma;
using namespace std;


SEXP mmIgroupDP(SEXP yvec, SEXP Xmat, SEXP Zmat, SEXP Wcase, SEXP Wperson, SEXP Mmat,
        SEXP personsvec, SEXP niterInt, SEXP nburnInt, SEXP nthinInt, SEXP ustrengthd,
        SEXP shapealph, SEXP ratebeta)
{
BEGIN_RCPP
    // Time run
    //wall_clock timer;
    //timer.tic();
    // Wrap SEXP objects in Rcpp Vector and Matrix objects.  Not copies,
    // but wrappers.
    // Data: observations, fe, re, mm weights, car weights, session groups,
    // client id's
    // Counters (number of); cbt cases, mcmc iter, burnin.
    NumericVector yr(yvec);
    NumericMatrix Xr(Xmat);
    NumericMatrix Zr(Zmat);
    NumericMatrix Wcr(Wcase);
    NumericMatrix Mr(Mmat);
    NumericMatrix Wp(Wperson);
    IntegerVector pr(personsvec);
    int niter = as<int>(niterInt);
    int nburn = as<int>(nburnInt);
    int nthin = as<int>(nthinInt);
    double ac = as<double>(shapealph);
    double bc = as<double>(ratebeta);
    double ustrength = as<double>(ustrengthd);
    int nkeep = (niter -  nburn)/nthin;

    // Extract key row and column dimensions we will need to sample parameters
    int nc =  Xr.nrow(), nf = Xr.ncol(), nr = Zr.ncol(), ns = Wcr.ncol();
    int npcbt = Wp.nrow(), ng = Mr.ncol();

    // Compute np, number of unique clients
    // algorithm assumes cases clustered by unique client
    IntegerVector dpr = diff(pr);
    Col<int32_t> diffpersons(dpr.begin(), nc - 1, false);
    int np = accu(diffpersons) + 1;
    /* int np = accumulate(dpr.begin,dpr.end,1); */

    // Create armadillo structures we will use for our computations
    // We set copy_aux_mem = false from the constructor, meaning we re-use
    // the memory to in which our object is stored
    // Xr.begin() is a pointer, presumably of type double*
    mat xmat(Xr.begin(), nc, nf, false);
    mat zmat(Zr.begin(), nc, nr, false);
    mat mmat(Mr.begin(), ns, ng, false);
    mat wcase(Wcr.begin(),nc, ns, false);
    mat wpers(Wp.begin(),npcbt, ns, false); /* mmwt */
    colvec y(yr.begin(), nc, false);
    Col<int32_t> persons(pr.begin(), nc, false);
    // set up data objects indexed by clustering unit
    field<mat> zsplit(np,1); field<mat> ztzsplit(np,1);
    field<mat> xsplit(np,1); field<mat> wsplit(np,1);
    field<colvec> ysplit(np,1); field<colvec> ytilsplit(np,1);
    CovUnitSplitMM(zmat, xmat, wcase, y, zsplit, ztzsplit, xsplit, wsplit, ysplit, persons);
    mat xtx; prodMatOne(xmat,xtx);  /* sampling fixed effects */
    mat wtwmat; prodMatOne(wcase,wtwmat); /* joint sampling of MM effects */
    mat mtmat; prodMatOne(mmat,mtmat); /* sampling of by-group MM effects mean */

    // Set random number generator state
    RNGScope scope; /* Rcpp */
    //arma_rng::set_seed_random(); /* arma */

    // Armadillo structures to capture results
    // Will be wrapped into a list element on return
    colvec Deviance(nkeep); Deviance.zeros();
    colvec Alpha(nkeep);
    colvec Conc(nkeep);
    mat Devmarg(nkeep,nc);
    mat Resid(nkeep,nc);
    mat Beta(nkeep,nf);
    mat B(nkeep,nr*np);
    rowvec brow(np*nr); brow.zeros();
    mat MM(nkeep,npcbt);
    mat Eta(nkeep,ng);
    mat U(nkeep,ns);
    mat Taub(nkeep,nr);
    colvec Tauu(nkeep);
    colvec Taue(nkeep);
    colvec Taueta(nkeep);
    icolvec numM(nkeep);
    imat S(nkeep,np);
    field<icolvec> Num(nkeep,1);
    field<ucolvec> bigS;
    field<mat> optPartition(3,1); /* Hold empirical probability of co-clustering matrix, index of L-sq clustering scores objects, and cluster identifiers per iteration */

    // Initialize parameter values
    /* cluster capture variables */
    double tauu, taubeta, taualph, taue, taueta;
    tauu = taubeta = taualph = taue = taueta = 1;
    rowvec taub(nr); taub.ones();
    IntegerVector ss = seq_len(np);
    icolvec s(ss.begin(), np, false);
    s -= 1; s = shuffle(s);
    int M = np; /* each person initialized with their own cluster */
    icolvec num(M); num.ones();
    mat    bstarmat = randn<mat>(M,nr)*sqrt(1/taub(0));
    mat    bmat(np,nr); bmat.zeros();
    colvec zb(nc); zb.zeros();
    double conc = 1; /* DP concentration parameter */
    ucolvec ordscore; ordscore.zeros(); /* LS score order */
    mat phat(np,np); phat.zeros(); /* empirical co-clustering probability matrix */
    /* remaining parameter set */
    colvec beta = randn<colvec>(nf)*sqrt(1/taubeta);
    colvec eta(ng); eta.zeros(); /* group means of session effects */
    colvec u = randn<colvec>(ns)*sqrt(1/tauu);
    colvec resid(nc); resid.zeros(); /* regression residual */
    colvec mm(npcbt); mm.zeros(); /* client effect mapped from session eff */
    mat pbmat(nr,nr); pbmat.zeros(); /* prior precision matrix for B */
    double alpha = rnorm( 1, 0, sqrt(1/taualph) )[0];
    /* goodness of fit related measures */
    double deviance; colvec devres(4); rowvec devmarg(nc);
    rowvec logcpo(nc); logcpo.zeros(); double lpml;
    /* capture samples to return */
    int oo = 0, kk;

    // set hyperparameter values
    double a1, a2, a4, a7, a8, b1, b2, b4, b7, b8;
    a1 = b1 = ustrength; /* tauu */
    a2 = b2 = 1; /* taub */
    a4 = b4 = 1; /* taue */
    a7 = b7 = 1; /* taueta*/
    a8 = ac; b8 = bc; /* c */

    // Conduct posterior samples and store results
    int k, l;
    for(k = 0; k < niter; k++)
    {
        //if( (k % 1000) == 0 ) cout << "Interation: " << k << endl;
        clustermmstep(zsplit, ztzsplit, xsplit, wsplit, ysplit, ytilsplit, pbmat, 
                u, beta, alpha,
                taue, taub, bstarmat, s, num, M, conc,
                nr);
        concstep(conc, M, np, a8, b8);
        bstarstep(zsplit, ztzsplit, ytilsplit, taue, pbmat, s, num, bstarmat,
                zb, bmat, np, nr);
        betalscholstep(xmat, xtx, wcase, y, beta, u, alpha, taue, zb, nf);
        uindetastep(xmat, wcase, wtwmat, wpers, mmat, beta, zb, y, eta, u, mm,
                alpha, taue, tauu, ns);
        etastep(mmat, mtmat, u, eta, tauu, taueta);
        alphalsstep(xmat, wcase, beta, zb, y, u, resid, alpha, taue, nc);
        taummigdpstep(bstarmat, mmat, resid, u, eta, tauu, taub, taue,
                taueta, a1, a2, a4, a7, b1, b2, b4, b7, M, ns, nc, ng);
        if(k >= nburn)
        {
            deviance =  dev(resid, taue); /* scalar double deviance */
            dmarg(resid, taue, devmarg); /* 1 x nc vector of densities*/
            kk = k - nburn;
            if( kk == ((oo+1)*nthin - 1) )
            {
                Deviance(oo) = deviance;
                Devmarg.row(oo)  = devmarg;
                Beta.row(oo) = trans(beta);
                Alpha(oo) = alpha;
                Conc(oo) = conc;
                for( l = 0; l < nr; l++)
                {
                    brow.cols( np*l,(np*(l+1)-1) ) = trans( bmat.col(l) );
                }
                B.row(oo) = brow;
                U.row(oo) = trans(u);
                MM.row(oo) = trans(mm);
                Tauu(oo) = tauu;
                Taue(oo) = taue;
                Taueta(oo) = taueta;
                Taub.row(oo) = taub;
                Resid.row(oo) = trans(resid);
                numM(oo) = M;
                S.row(oo) = trans(s);
                Num(oo,0) = num;
                oo += 1; /* increment return sample counter */
            } /* end conditional statement returning samples */
        } /* end if k > burnin, record results for return */

    } /* end MCMC loop over k */

    // compute least squares cluster
    lsqcluster(S, Num, ordscore, phat, bigS);
    optPartition(0,0) = phat;
    optPartition(1,0) = conv_to<mat>::from(ordscore);
    optPartition(2,0) = conv_to<mat>::from(S);
    // DIC
    dic3comp(Deviance, Devmarg, devres); /* devres = c(dic,dbar,dhat,pd) */
    cpo(Devmarg, logcpo, lpml);

    // Return results
    return Rcpp::List::create(Rcpp::Named("Deviance") = Deviance,
                                  Rcpp::Named("devres") = devres,
                                  //Rcpp::Named("logcpo") = logcpo,
                                  Rcpp::Named("lpml") = lpml,
				  Rcpp::Named("Beta") = Beta,
				  Rcpp::Named("Alpha") = Alpha,
                                  Rcpp::Named("C") = Conc,
                                  Rcpp::Named("B") = B,
                                  Rcpp::Named("U") = U,
                                  Rcpp::Named("Eta") = Eta,
                                  Rcpp::Named("MM") = MM,
                                  Rcpp::Named("Tauu") = Tauu,
                                  Rcpp::Named("Taue") = Taue,
                                  //Rcpp::Named("Taueta") = Taueta,
                                  Rcpp::Named("Taub") = Taub,
                                  Rcpp::Named("Residuals") = Resid,
                                  Rcpp::Named("M") = numM,
                                  //Rcpp::Named("S") = S,
                                  Rcpp::Named("Num") = Num,
                                  Rcpp::Named("optPartition") = optPartition,
                                  Rcpp::Named("bigSmin") = bigS
                                  );

    // Print CPU runtime
     //double n_secs = timer.toc();
     //cout << "took " << n_secs << " seconds" << endl;
END_RCPP
} /* end MCMC function returning SEXP */

SEXP taummigdpstep(const mat& bstarmat, const mat& mmat, const colvec& resid,
            const colvec& u, const colvec& eta, double& tauu, rowvec& taub, 
            double& taue, double& taueta, double a1, double a2, double a4, double a7,
            double b1, double b2, double b4, double b7, int M, int ns, int nc, int ng)
    {
        BEGIN_RCPP
        double a, b;
        /* tauu */
        a = 0.5*double(ns) + a1;
        colvec uresid = u - mmat*eta;
        b = 0.5* dot( uresid , uresid ) + b1;
        tauu = rgamma(1, a, (1/b))[0];
        /* tau0 */
        int nr = bstarmat.n_cols;
        a = 0.5*double(M) + a2;
        int i;
        for(i = 0; i < nr; i++)
        {
            b =  0.5*dot( bstarmat.col(i), bstarmat.col(i) ) + b2;
            taub(i) = rgamma(1, a, (1/b))[0];
        }
        /* taue */
        a = 0.5*double(nc) + a4;
        b = 0.5*dot( resid, resid ) + b4;
        taue = rgamma(1, a, (1/b))[0];
        /* taueta */
        a = 0.5*double(ng) + a7;
        b = 0.5* dot( eta , eta ) + b7;
        taueta = rgamma(1, a, (1/b))[0];
        END_RCPP
    } /* end function taustep to sample precision parameters */
