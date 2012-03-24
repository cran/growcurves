#include "mmC.h"

using namespace Rcpp;
using namespace arma;
using namespace std;


SEXP mmIplusDP(SEXP yvec, SEXP Xmat, SEXP Zmat, SEXP Wcase, SEXP Wperson,
        SEXP personsvec, SEXP niterInt, SEXP nburnInt, SEXP nthinInt, 
        SEXP ustrengthd, SEXP shapealph)
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
    NumericMatrix Wp(Wperson);
    IntegerVector pr(personsvec);
    int niter = as<int>(niterInt);
    int nburn = as<int>(nburnInt);
    int nthin = as<int>(nthinInt);
    double ac = as<double>(shapealph);
    double ustrength = as<double>(ustrengthd);
    int nkeep = (niter -  nburn)/nthin;

    // Extract key row and column dimensions we will need to sample parameters
    int nc =  Xr.nrow(), nf = Xr.ncol(), nr = Zr.ncol(), ns = Wcr.ncol();
    int npcbt = Wp.nrow();

    // Compute np, number of unique clients
    // algorithm assumes cases clustered by unique client
    IntegerVector dpr = diff(pr);
    icolvec diffpersons(dpr.begin(), nc - 1, false);
    int np = accu(diffpersons) + 1;
    /* int np = accumulate(dpr.begin,dpr.end,1); */

    // Create armadillo structures we will use for our computations
    // We set copy_aux_mem = false from the constructor, meaning we re-use
    // the memory to in which our object is stored
    // Xr.begin() is a pointer, presumably of type double*
    mat xmat(Xr.begin(), nc, nf, false);
    mat zmat(Zr.begin(), nc, nr, false);
    mat wcase(Wcr.begin(),nc, ns, false);
    mat wpers(Wp.begin(),npcbt, ns, false); /* mmwt */
    colvec y(yr.begin(), nc, false);
    icolvec persons(pr.begin(), nc, false);
    field<mat> zsplit(np,1);
    field<colvec> ytilsplit(np,1);

    // Set random number generator state
    RNGScope scope; /* Rcpp */
    srand ( time(NULL) ); /* arma */

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
    mat U(nkeep,ns);
    mat Taub(nkeep,nr);
    colvec Tauu(nkeep);
    colvec Taue(nkeep);
    icolvec numM(nkeep);
    imat S(nkeep,np);
    field<icolvec> Num(nkeep,1);
    field<ucolvec> bigS;

    // Initialize parameter values
    /* cluster capture variables */
    double tauu, taubeta, taualph, taue;
    tauu = taubeta = taualph = taue = 1;
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
    /* remaining parameter set */
    colvec beta = randn<colvec>(nf)*sqrt(1/taubeta);
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
    double a1, a2, a4, a7, b1, b2, b4, b7;
    a1 = b1 = ustrength; /* tauu */
    a2 = b2 = 1; /* taub */
    a4 = b4 = 1; /* taue */
    a7 = ac; b7 = 1; /* c */

    // Conduct posterior samples and store results
    int k, l;
    for(k = 0; k < niter; k++)
    {
        //if( (k % 1000) == 0 ) cout << "Interation: " << k << endl;
        clustermmstep(xmat, zmat, zsplit, ytilsplit, pbmat, y,
                wcase, u, beta, alpha,
                taue, taub, persons, bstarmat, s, num, M, conc,
                np, nr);
        concstep(conc, M, np, a7, b7);
        bstarstep(zsplit, ytilsplit, taue, pbmat, s, num, bstarmat,
                zb, bmat, np, nr);
        betalscholstep(xmat, wcase, y, beta, u, alpha, taue, zb, nf);
        uindstep(xmat, wcase, wpers, beta, zb, y, u, mm,
                alpha, taue, tauu, ns);
        alphalsstep(xmat, wcase, beta, zb, y, u, resid, alpha, taue, nc);
        taummidpstep(bstarmat, resid, u, tauu, taub, taue, a1, a2, a4, b1, b2,
                b4, M, ns, nc);
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
                Taub.row(oo) = taub;
                Resid.row(oo) = trans(resid);
                numM(oo) = M;
                S.row(oo) = trans(s);
                Num(oo,0) = num;
                oo += 1; /*increment sample return counter */
            } /* end conditional statement on returning samples */
        } /* end if k > burnin, record results for return */

    } /* end MCMC loop over k */

    // compute least squares cluster
    lsqcluster(S, Num, ordscore, bigS);
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
                                  Rcpp::Named("MM") = MM,
                                  Rcpp::Named("Tauu") = Tauu,
                                  Rcpp::Named("Taue") = Taue,
                                  Rcpp::Named("Taub") = Taub,
                                  Rcpp::Named("Residuals") = Resid,
                                  Rcpp::Named("M") = numM,
                                  Rcpp::Named("S") = S,
                                  Rcpp::Named("Num") = Num,
                                  //Rcpp::Named("ordscore") = ordscore,
                                  Rcpp::Named("bigSmin") = bigS
                                  );

    // Print CPU runtime
     //double n_secs = timer.toc();
     //cout << "took " << n_secs << " seconds" << endl;
END_RCPP
} /* end MCMC function returning SEXP */


    SEXP taummidpstep(const mat& bstarmat, const colvec& resid,
            const colvec& u, double& tauu, rowvec& taub, double& taue,
            double a1, double a2, double a4, double b1, double b2, double b4,
            int M, int ns, int nc)
    {
        BEGIN_RCPP
        double a, b;
        /* tauu */
        a = 0.5*double(ns) + a1;
        b = 0.5* dot( u , u ) + b1;
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
        END_RCPP
    } /* end function taustep to sample precision parameters */
